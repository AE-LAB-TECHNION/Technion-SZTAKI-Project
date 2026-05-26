# AGENTS.md

## Principles
- Prefer minimal changes
- Avoid unnecessary abstractions
- Do not rewrite working numerical code
- Preserve existing function signatures
- Keep scripts readable and vectorized where appropriate
- Do not introduce classes unless explicitly requested

## Workflow
- Explain intended changes before editing
- Modify only files directly related to the task
- Avoid touching plotting/visualization code unless necessary
- Prefer incremental edits over large refactors

## Validation
- Ensure MATLAB syntax is valid
- Check for dimension mismatches
- Preserve current outputs unless explicitly changing behavior

## MATLAB Style
- Use clear variable names
- Avoid deeply nested logic
- Prefer functions over large monolithic scripts
- Keep compatibility with current MATLAB version
