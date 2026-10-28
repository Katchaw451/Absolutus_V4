#!/usr/bin/env bash
# Initialise git, commit, create the PUBLIC GitHub repo and push. Uses your existing git/gh credentials;
# it never reads, prints or stores tokens.
set -euo pipefail
REPO="${1:-Katchaw451/Absolutus_V4}"
cd "$(dirname "$0")/.."
command -v git >/dev/null || { echo "git not found"; exit 1; }
git config user.name  >/dev/null || { echo "Set git user.name first:  git config --global user.name 'Your Name'"; exit 1; }
git config user.email >/dev/null || { echo "Set git user.email first: git config --global user.email you@example.com"; exit 1; }
[ -d .git ] || git init -b main
git add -A
git diff --cached --quiet || git commit -m "Absolutus V4 v0.1.0: validation-stage pipeline (simulated data only)"
if git remote get-url origin >/dev/null 2>&1; then
  git push -u origin main
elif command -v gh >/dev/null && gh auth status >/dev/null 2>&1; then
  gh repo create "$REPO" --public --source=. --remote=origin --push
else
  echo "gh CLI missing or not authenticated. Create an EMPTY public repo named $REPO on GitHub, then run:"
  echo "  git remote add origin git@github.com:$REPO.git   # or the https URL"
  echo "  git push -u origin main"
  exit 2
fi
