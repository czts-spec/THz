# .chatman/README

This project is initialized for `chatman` usage and follows the profile-level rules from `czts-spec/dtba3a-del`.

The project-local instructions live in `.chatman/rules.md`; this file is the pointer for the repository itself, so the project remains aware of the tool without copying the tool sources into the repo.

## Maintained conventions

- Only one branch is used: `main`.
- Project-specific rules stay in `.chatman/rules.md`.
- Chat logs are managed through `chatman` and remain outside the repository unless explicitly exported.
- The repository follows the profile ban on `Co-Authored-By`, tool signatures, and claims of authorship.

## Commands

```bash
chatman rules
chatman repos
chatman find --prompts --grep '...'
chatman search '...'
chatman export
```

The `chatman` tool is the canonical interface for the project log and rules context.
