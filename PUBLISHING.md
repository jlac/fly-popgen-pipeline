# Publish this repository to GitHub

Target repository: **public** `jlac/fly-popgen-pipeline`.

The source archive excludes original sequencing files, real ICP measurements, and individual-level metadata. Review its contents before posting anything publicly; public GitHub repositories are visible to everyone.

## Option A: GitHub CLI (recommended)

Install [GitHub CLI](https://cli.github.com/), sign in to your `jlac` account, and from the extracted `fly-popgen-pipeline` folder run:

```bash
gh auth login
git init -b main
git add .
git status --short
git commit -m "Initial public release of Fly PopGen Pipeline"
gh repo create jlac/fly-popgen-pipeline --public --source=. --remote=origin --push
```

If Git asks for an author identity, set it with `git config user.name` and `git config user.email`, then repeat the commit.

## Option B: GitHub website

1. Visit https://github.com/new while logged in as `jlac`.
2. Choose repository name `fly-popgen-pipeline`, select **Public**, and create it **without** initializing a README, license, or `.gitignore`.
3. From this extracted folder run:

```bash
git init -b main
git add .
git status --short
git commit -m "Initial public release of Fly PopGen Pipeline"
git remote add origin https://github.com/jlac/fly-popgen-pipeline.git
git push -u origin main
```

GitHub's UI and CLI may prompt for authentication and permissions.

## After publishing

- Check that `README.md`, `fly_fst_pipeline/`, `tools/`, and `.github/workflows/tests.yml` are present.
- Check the Actions tab for results; a synthetic test currently skips under some environments.
- Decide whether to add a software license; this distribution intentionally does not presume one.
- Keep real genetic and phenotype data outside the public repo.
