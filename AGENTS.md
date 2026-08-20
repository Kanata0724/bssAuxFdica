# AGENTS.md

## Project purpose
This repository is used for research on blind source separation and AuxFDICA-related experiments.

## Working rules
- Before starting work, read `PROGRESS.md` and understand the current research status.
- Do not modify files unrelated to the requested task.
- Do not delete existing experimental code or results without explicit permission.
- Prefer small, reviewable changes.
- Explain significant changes before or after making them.

## Experiments
- Experimental outputs should be saved under `output/`.
- Do not add files under `output/` to Git.
- Record important experimental conditions and results in `PROGRESS.md`.
- When possible, record:
  - purpose of the experiment
  - changed parameters
  - executed script or command
  - important numerical results
  - interpretation
  - next action

## Git
- Do not push to GitHub unless explicitly instructed.
- Do not force-push.
- Do not rewrite Git history.
- Before committing, check `git status` and confirm that no large or unintended files are included.
- Never commit secrets, credentials, API keys, or private data.

## Safety
- Do not modify or overwrite datasets unless explicitly instructed.
- Ask before performing destructive operations.