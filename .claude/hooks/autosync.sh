#!/usr/bin/env bash
# Stop hook: commit and push any working-tree changes after each Claude turn.
# Never blocks Claude: all failures are logged to .claude/hooks/autosync.log and exit 0.
cd "$CLAUDE_PROJECT_DIR" || exit 0
LOG=.claude/hooks/autosync.log

# Skip while a merge/rebase is in progress or if Dropbox left lock files.
if [ -d .git/rebase-merge ] || [ -d .git/rebase-apply ] || [ -f .git/MERGE_HEAD ] || [ -f .git/index.lock ]; then
  echo "$(date '+%F %T') skipped: git busy" >> "$LOG"; exit 0
fi

git add -A -- . ':!*.DS_Store' ':!.claude/hooks/autosync.log' >> "$LOG" 2>&1
if git diff --cached --quiet; then exit 0; fi

git commit -q -m "Auto-sync: $(git diff --cached --name-only | head -5 | paste -sd, -)" \
  -m "Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>" >> "$LOG" 2>&1 \
  || { echo "$(date '+%F %T') commit failed" >> "$LOG"; exit 0; }

git pull --rebase -q >> "$LOG" 2>&1 || { git rebase --abort >> "$LOG" 2>&1; echo "$(date '+%F %T') pull --rebase failed, not pushed" >> "$LOG"; exit 0; }
git push -q >> "$LOG" 2>&1 || echo "$(date '+%F %T') push failed" >> "$LOG"
exit 0
