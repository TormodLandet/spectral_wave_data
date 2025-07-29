import subprocess
import os


def should_run_quick_tests():
    """
    Determine if quick tests should be run based on the environment variable.
    """
    return os.environ.get("SWD_TEST_TYPE", "normal") == "quick"


def run_swd_meta_check(file_swd, **tags):
    """
    Run the swd_meta command and check the output for specific tags
    that should be present in the SWD file metadata.
    """
    result = subprocess.run(["swd_meta", file_swd], capture_output=True)

    if result.returncode != 0:
        print(f"swd_meta failed with return code {result.returncode}")
        print(f"stdout:\n{'-' * 80}\n{result.stdout.decode()}\n{'-' * 80}")
        print(f"stderr:\n{'-' * 80}\n{result.stderr.decode()}\n{'-' * 80}")
        raise Exception(f"swd_meta failed with return code {result.returncode}")

    text = result.stdout.decode()

    def check(tag, val):
        text_ok = f"{tag + ':':<8} {val}"
        assert text_ok in text, f"Missing text {text_ok}, got:\n{'-' * 80}\n{text}"

    for tag, val in tags.items():
        check(tag, val)
