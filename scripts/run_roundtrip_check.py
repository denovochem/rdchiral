"""Run roundtrip consistency check across different rdchiral environments."""

import argparse
import os
import random
import shutil
import subprocess
import tempfile
from pathlib import Path


def _run(
    cmd: list[str], *, env: dict[str, str] | None = None, cwd: Path | None = None
) -> None:
    printable = " ".join(cmd)
    print(f"\n$ {printable}")
    subprocess.run(cmd, check=True, env=env, cwd=str(cwd) if cwd is not None else None)


def _venv_python(venv_dir: Path) -> Path:
    if os.name == "nt":
        return venv_dir / "Scripts" / "python.exe"
    return venv_dir / "bin" / "python"


def _build_env_from_url(*, install_spec: str, venv_dir: Path, reinstall: bool) -> None:
    """Build a uv venv and install a package from a pip-compatible spec."""
    if reinstall and venv_dir.exists():
        shutil.rmtree(venv_dir)

    if not venv_dir.exists():
        _run(["uv", "venv", str(venv_dir)])

    venv_python = _venv_python(venv_dir)
    _run(
        [
            "uv",
            "pip",
            "install",
            "--python",
            str(venv_python),
            install_spec,
            "rdkit",
        ],
    )


def _conda_python(env_dir: Path) -> Path:
    if os.name == "nt":
        return env_dir / "python.exe"
    return env_dir / "bin" / "python"


def _find_conda() -> str:
    """Resolve the full path to conda/mamba/micromamba."""
    for name in ("conda", "mamba", "micromamba"):
        found = shutil.which(name)
        if found:
            return found

    home = Path.home()
    for candidate in (
        home / "miniforge3" / "bin" / "conda",
        home / "mambaforge" / "bin" / "conda",
        home / "miniconda3" / "bin" / "conda",
        home / "anaconda3" / "bin" / "conda",
    ):
        if candidate.exists():
            return str(candidate)

    try:
        result = subprocess.run(
            ["bash", "-i", "-c", "which conda"],
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass

    raise FileNotFoundError(
        "Could not find conda, mamba, or micromamba. "
        "Please install one of them or ensure it is on PATH."
    )


def _build_conda_env(*, env_dir: Path, reinstall: bool) -> None:
    """Create a conda prefix env and install rdchiral_cpp from conda-forge."""
    conda = _find_conda()

    if reinstall and env_dir.exists():
        shutil.rmtree(env_dir)

    if not env_dir.exists():
        _run(
            [
                conda,
                "create",
                "--prefix",
                str(env_dir),
                "-c",
                "conda-forge",
                "rdchiral_cpp",
                "-y",
            ]
        )
    else:
        _run(
            [
                conda,
                "install",
                "--prefix",
                str(env_dir),
                "-c",
                "conda-forge",
                "rdchiral_cpp",
                "-y",
            ]
        )


def _build_env(
    *, repo_root: Path, venv_dir: Path, use_mypyc: bool, reinstall: bool
) -> None:
    if reinstall and venv_dir.exists():
        shutil.rmtree(venv_dir)

    if not venv_dir.exists():
        _run(["uv", "venv", str(venv_dir)])

    venv_python = _venv_python(venv_dir)

    env = os.environ.copy()
    env["RDCHIRAL_USE_MYPYC"] = "1" if use_mypyc else "0"

    _run(
        ["uv", "pip", "install", "--python", str(venv_python), "."],
        env=env,
        cwd=repo_root,
    )

    _run(
        [
            "uv",
            "pip",
            "install",
            "--python",
            str(venv_python),
            "rdkit",
        ],
        env=env,
    )


def _run_roundtrip_check(
    *,
    python: Path,
    repo_root: Path,
    check_script_path: Path,
) -> int:
    """Run the roundtrip check and return the consistent count."""
    # Critical: run from a directory that does NOT contain the repo to avoid importing
    # the in-tree rdchiral/*.py instead of the installed package.
    with tempfile.TemporaryDirectory(prefix="rdchiral_roundtrip_") as d:
        tmpdir = Path(d)
        local_script = tmpdir / check_script_path.name
        shutil.copy2(check_script_path, local_script)
        cmd = [str(python), str(local_script)]

        # Pass repository root as environment variable
        env = os.environ.copy()
        env["RDCHIRAL_REPO_ROOT"] = str(repo_root) + "/scripts"

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=tmpdir,
            env=env,
        )
        if result.returncode != 0:
            print(f"Error running check script: {result.stderr}")
            raise subprocess.CalledProcessError(
                result.returncode, cmd, output=result.stdout, stderr=result.stderr
            )
        return int(result.stdout.strip())


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run roundtrip consistency check across rdchiral environments"
    )
    parser.add_argument(
        "--venv-py",
        default=".venv-py",
        help="Path for the pure-Python venv (default: .venv-py)",
    )
    parser.add_argument(
        "--venv-mypyc",
        default=".venv-mypyc",
        help="Path for the mypyc venv (default: .venv-mypyc)",
    )
    parser.add_argument(
        "--venv-default",
        default=".venv-default",
        help="Path for the original rdchiral venv (default: .venv-default)",
    )
    parser.add_argument(
        "--venv-cpp",
        default=".conda-rdchiral-cpp",
        help="Path for the rdchiral_cpp conda env (default: .conda-rdchiral-cpp)",
    )
    parser.add_argument(
        "--reinstall",
        action="store_true",
        help="Delete and recreate venvs before installing",
    )
    parser.add_argument(
        "--shuffle-seed",
        type=int,
        default=None,
        help="Randomize the order of environment checks",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parent.parent
    check_script_path = (repo_root / "scripts" / "roundtrip_check_script.py").resolve()
    if not check_script_path.exists():
        raise FileNotFoundError(f"Check script not found: {check_script_path}")

    venv_py = (repo_root / args.venv_py).resolve()
    venv_mypyc = (repo_root / args.venv_mypyc).resolve()
    venv_default = (repo_root / args.venv_default).resolve()
    venv_cpp = (repo_root / args.venv_cpp).resolve()

    results: dict[str, int] = {}

    def _check_pure_python() -> None:
        print("\n=== Pure Python ===")
        _build_env(
            repo_root=repo_root,
            venv_dir=venv_py,
            use_mypyc=False,
            reinstall=args.reinstall,
        )
        py_python = _venv_python(venv_py)
        consistent = _run_roundtrip_check(
            python=py_python,
            repo_root=repo_root,
            check_script_path=check_script_path,
        )
        results["pure_python"] = consistent
        print(f"Consistent: {consistent}")

    def _check_mypyc() -> None:
        print("\n=== MYPYC ===")
        _build_env(
            repo_root=repo_root,
            venv_dir=venv_mypyc,
            use_mypyc=True,
            reinstall=args.reinstall,
        )
        mypyc_python = _venv_python(venv_mypyc)
        consistent = _run_roundtrip_check(
            python=mypyc_python,
            repo_root=repo_root,
            check_script_path=check_script_path,
        )
        results["mypyc"] = consistent
        print(f"Consistent: {consistent}")

    def _check_default() -> None:
        print("\n=== Default (pip) ===")
        _build_env_from_url(
            install_spec="rdchiral",
            venv_dir=venv_default,
            reinstall=args.reinstall,
        )
        default_python = _venv_python(venv_default)
        consistent = _run_roundtrip_check(
            python=default_python,
            repo_root=repo_root,
            check_script_path=check_script_path,
        )
        results["default"] = consistent
        print(f"Consistent: {consistent}")

    def _check_cpp() -> None:
        print("\n=== C++ (rdchiral_cpp) ===")
        _build_conda_env(
            env_dir=venv_cpp,
            reinstall=args.reinstall,
        )
        cpp_python = _conda_python(venv_cpp)
        consistent = _run_roundtrip_check(
            python=cpp_python,
            repo_root=repo_root,
            check_script_path=check_script_path,
        )
        results["cpp"] = consistent
        print(f"Consistent: {consistent}")

    checks = [
        ("pure_python", _check_pure_python),
        ("mypyc", _check_mypyc),
        ("default", _check_default),
        ("cpp", _check_cpp),
    ]

    if args.shuffle_seed is not None:
        random.seed(args.shuffle_seed)
        random.shuffle(checks)

    for name, check_fn in checks:
        check_fn()

    print("\n=== Summary ===")
    for name, consistent in results.items():
        print(f"  {name}: {consistent}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
