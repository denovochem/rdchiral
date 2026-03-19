import argparse
import os
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

    # `uv` is typically installed globally, not inside the venv.
    # Use `--python` to ensure the install targets this venv.
    _run(
        ["uv", "pip", "install", "--python", str(venv_python), ".", "-v"],
        env=env,
        cwd=repo_root,
    )


def _verify_import(*, python: Path) -> None:
    with tempfile.TemporaryDirectory(prefix="rdchiral_importcheck_") as d:
        tmpdir = Path(d)
        _run(
            [
                str(python),
                "-c",
                "import rdchiral.main; print(rdchiral.main.__file__)",
            ],
            cwd=tmpdir,
        )


def _run_benchmark(*, python: Path, repo_root: Path, benchmark_path: Path) -> None:
    # Critical: run from a directory that does NOT contain the repo to avoid importing
    # the in-tree rdchiral/*.py instead of the installed package.
    with tempfile.TemporaryDirectory(prefix="rdchiral_bench_") as d:
        tmpdir = Path(d)
        local_benchmark = tmpdir / benchmark_path.name
        shutil.copy2(benchmark_path, local_benchmark)
        _run([str(python), str(local_benchmark)], cwd=tmpdir)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--benchmark",
        default="speed_benchmark_script.py",
        help="Path to the benchmark script (default: speed_benchmark_script.py)",
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
        "--reinstall",
        action="store_true",
        help="Delete and recreate venvs before installing",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parent
    benchmark_path = (repo_root / args.benchmark).resolve()
    if not benchmark_path.exists():
        raise FileNotFoundError(f"Benchmark script not found: {benchmark_path}")

    venv_py = (repo_root / args.venv_py).resolve()
    venv_mypyc = (repo_root / args.venv_mypyc).resolve()

    print("=== Building pure-Python environment ===")
    _build_env(
        repo_root=repo_root, venv_dir=venv_py, use_mypyc=False, reinstall=args.reinstall
    )
    py_python = _venv_python(venv_py)
    print("--- Import verification (pure python) ---")
    _verify_import(python=py_python)
    print("--- Running benchmark (pure python) ---")
    _run_benchmark(python=py_python, repo_root=repo_root, benchmark_path=benchmark_path)

    print("\n=== Building mypyc environment ===")
    _build_env(
        repo_root=repo_root,
        venv_dir=venv_mypyc,
        use_mypyc=True,
        reinstall=args.reinstall,
    )
    mypyc_python = _venv_python(venv_mypyc)
    print("--- Import verification (mypyc) ---")
    _verify_import(python=mypyc_python)
    print("--- Running benchmark (mypyc) ---")
    _run_benchmark(
        python=mypyc_python, repo_root=repo_root, benchmark_path=benchmark_path
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
