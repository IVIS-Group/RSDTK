"""
Pre-push audit for RSDTK.

Scans the repository for things that should not be made public, or that will
bloat the repository permanently once committed. Run this from the repository
root before the first push, and ideally before any release.

    python scripts/pre_push_audit.py

Exits 0 if nothing needs attention, 1 if there are findings to review.
Findings are advisory, not definitive — read them and judge.
"""

import os
import re
import sys

# ---------------------------------------------------------------------------
# What to look for
# ---------------------------------------------------------------------------

# Absolute paths. Requires a backslash separator so that URL schemes such as
# "https://" are not matched. The CLI help text in the converters legitimately
# contains generic examples like D:\Data\AIRS, reported separately as info.
ABS_PATH = re.compile(r"(?<![\w/])[A-Za-z]:\\{1,2}[\w\s.\\-]+")
GENERIC_EXAMPLE = re.compile(r"[A-Za-z]:\\{1,2}Data\\{1,2}\w+")

UNIX_ABS_PATH = re.compile(r"[\"'](/(?:home|Users|mnt|media)/[\w./-]+)[\"']")

CREDENTIAL_PATTERNS = [
    (re.compile(r"(?i)\bpass(?:word|wd)\s*=\s*[\"'][^\"']+[\"']"), "password assignment"),
    (re.compile(r"(?i)\bapi[_-]?key\s*=\s*[\"'][^\"']+[\"']"), "API key assignment"),
    (re.compile(r"(?i)\b(?:auth|access|bearer)?_?token\s*=\s*[\"'][^\"']+[\"']"), "token assignment"),
    (re.compile(r"(?i)\bsecret\s*=\s*[\"'][^\"']+[\"']"), "secret assignment"),
    (re.compile(r"(?i)netrc"), "netrc reference"),
    (re.compile(r"(?i)urs\.earthdata\.nasa\.gov.*[\"']"), "Earthdata credential URL"),
    (re.compile(r"eyJ[A-Za-z0-9_-]{20,}"), "possible JWT"),
]

# Personal identifiers worth a look before going public
PERSONAL_PATTERNS = [
    (re.compile(r"(?i)\b[\w.+-]+@(?!example\.|test\.)[\w-]+\.[\w.]+"), "email address"),
    (re.compile(r"(?i)C:\\Users\\[^\\\s\"']+"), "Windows user directory"),
]

# File types that should never be committed
BANNED_EXTENSIONS = {
    ".hdf", ".he5", ".h5", ".nc", ".jp2", ".img",
    ".dat", ".hdr", ".met", ".exe", ".msi",
}
# .tif is banned except under docs/
CONDITIONAL_EXTENSIONS = {".tif", ".tiff"}

LARGE_FILE_MB = 5

SKIP_DIRS = {
    ".git", "__pycache__", ".pytest_cache", "venv", ".venv", "env",
    "node_modules", "build", "dist", ".idea", ".vscode",
}

TEXT_EXTENSIONS = {
    ".py", ".md", ".txt", ".yml", ".yaml", ".json", ".cfg", ".ini",
    ".toml", ".bat", ".sh", ".iss", ".html", ".css", ".js", ".cff",
}


# ---------------------------------------------------------------------------

class Findings:
    def __init__(self):
        self.critical = []
        self.review = []
        self.info = []

    def add_critical(self, path, line_no, message, excerpt=""):
        self.critical.append((path, line_no, message, excerpt))

    def add_review(self, path, line_no, message, excerpt=""):
        self.review.append((path, line_no, message, excerpt))

    def add_info(self, path, line_no, message, excerpt=""):
        self.info.append((path, line_no, message, excerpt))


def walk_files(root):
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
        for fn in filenames:
            yield os.path.join(dirpath, fn)


def check_file_types(path, root, findings):
    rel = os.path.relpath(path, root)
    ext = os.path.splitext(path)[1].lower()

    if ext in BANNED_EXTENSIONS:
        findings.add_critical(rel, None,
                              f"Data/binary file ({ext}) should not be committed")
        return

    if ext in CONDITIONAL_EXTENSIONS:
        if not rel.replace("\\", "/").startswith("docs/"):
            findings.add_critical(rel, None,
                                  f"Raster file ({ext}) outside docs/ should not be committed")
        return

    try:
        size_mb = os.path.getsize(path) / (1024 * 1024)
    except OSError:
        return
    if size_mb > LARGE_FILE_MB:
        findings.add_review(rel, None,
                            f"Large file: {size_mb:.1f} MB. Attach to a Release instead?")


# This script necessarily contains the patterns it searches for, so scanning it
# produces only self-matches.
SELF_EXCLUDE = {"pre_push_audit.py"}


def scan_text(path, root, findings):
    rel = os.path.relpath(path, root)
    ext = os.path.splitext(path)[1].lower()
    if ext not in TEXT_EXTENSIONS:
        return
    if os.path.basename(path) in SELF_EXCLUDE:
        return

    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            lines = fh.readlines()
    except OSError:
        return

    for i, line in enumerate(lines, 1):
        stripped = line.strip()

        for pattern, label in CREDENTIAL_PATTERNS:
            if pattern.search(line):
                findings.add_critical(rel, i, f"Possible {label}", stripped[:100])

        for pattern, label in PERSONAL_PATTERNS:
            m = pattern.search(line)
            if m:
                findings.add_review(rel, i, f"Possible {label}", stripped[:100])

        for m in ABS_PATH.finditer(line):
            excerpt = stripped[:100]
            if GENERIC_EXAMPLE.search(line):
                findings.add_info(rel, i, "Generic path example (probably fine)", excerpt)
            else:
                findings.add_review(rel, i, "Absolute Windows path", excerpt)

        m = UNIX_ABS_PATH.search(line)
        if m:
            findings.add_review(rel, i, "Absolute Unix path", stripped[:100])


def check_required_files(root, findings):
    required = [
        "README.md", "LICENSE", ".gitignore",
        "CITATION.cff", "requirements.txt",
    ]
    for fn in required:
        if not os.path.exists(os.path.join(root, fn)):
            findings.add_review(fn, None, "Expected file is missing")


def report(findings):
    def emit(items, heading, symbol):
        if not items:
            return
        print(f"\n{symbol} {heading} ({len(items)})")
        print("-" * 70)
        for path, line_no, message, excerpt in items:
            loc = f"{path}:{line_no}" if line_no else path
            print(f"  {loc}")
            print(f"      {message}")
            if excerpt:
                print(f"      > {excerpt}")

    emit(findings.critical, "MUST FIX before pushing", "[!]")
    emit(findings.review, "REVIEW these", "[?]")
    emit(findings.info, "Informational", "[i]")


def main():
    root = os.getcwd()
    print(f"Auditing: {root}")

    findings = Findings()

    for path in walk_files(root):
        check_file_types(path, root, findings)
        scan_text(path, root, findings)

    check_required_files(root, findings)

    report(findings)

    print("\n" + "=" * 70)
    n_crit = len(findings.critical)
    n_rev = len(findings.review)

    if n_crit == 0 and n_rev == 0:
        print("Nothing flagged. Note this is a pattern scan, not a guarantee.")
        return 0

    print(f"{n_crit} to fix, {n_rev} to review.")
    print()
    print("Reminder: anything already committed stays in git history even after")
    print("deletion. If a secret has been pushed, rotate it rather than relying")
    print("on removing the file.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
