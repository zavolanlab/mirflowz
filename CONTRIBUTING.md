# Guidelines for contributing

We would kindly ask you to abide by our [Code of Conduct][coc] in all
interactions with the community when contributing to this project, regardless
of the type of contribution. We will not accept any offensive or demeaning
behavior towards others and will take any necessary steps to ensure that
everyone is treated with respect and dignity.

## Submitting Issues

Please, use the [issue tracker][issue-tracker] to:

- find issues to work on
- report bugs
- propose features
- discuss future directions

Please, choose a template when submitting an issue: choose the
[**bug report** template][bug-report] only when reporting bugs; for all other
issues, choose the [**feature request** template][feat-request]. Please follow
the instructions in the templates.

You do not need to worry about adding labels or milestones for an issue, the
project maintainers will do that for you. However, it is important that all
issues are written concisely, yet with enough detail and with proper
references (links, screenshots, etc) to allow other contributors to start
working on them. For bug reports, it is essential that they include all
information required to reproduce the bug.

Please **do not** use the issue tracker to ask usage questions, installation
problems etc., unless they appear to be bugs. For these issues, please use
the [communication channels](#other-communication) outlined below.

## Other Communication

Send us an [email][contact] if you want to reach out to us regarding sensitive
matters such as Code of Conduct issues, personal inquiries or collaboration
requests.

For anything else, please use the [GitHub Discussions forum][discussion-forum].

## Contributing Code Changes

We are using [Git][git], [GitHub][github] and [Git Flow][git-flow].

> **Note:** If you are a **beginner** and do not have a lot of experience with
> this sort of workflow, please do not feel overwhelmed. We will guide you
> through the process until you feel comfortable using it. And do not worry
> about mistakes either - everybody does them. Often! Our project layout makes
> it very very hard for anyone to cause irreversible harm, so relax, try things
> out, take your time and enjoy the work! :)

### Code Style and Testing

To make it easier for everyone to maintain, read and contribute to the code,
as well as to ensure that the code base is robust and of high quality, we
would kindly ask you to stick to the following guidelines for code style and
testing.

- Please use a recent version of [Python 3][py] (3.12.0+)
- Please try to conform to the used code, docstrings and commenting style
  within a project to maintain consistency
- Please use [type hints][py-typing] for all function/method signatures
  (exception: tests)
- Please use the following linters (see configuration files in the repository's
  root directory, _e.g._ `pylint.cfg`, for settings):
  - [`flake8`][py-flake8]
  - [`mypy`][py-mypy] to help with type hints
  - [`pylint`][py-pylint] (use available [configuration][py-pylint-conf])
- Please use the following test suits:
  - [`pytest`][py-pytest]
  - [`coverage`][py-coverage]

### Commit Messages

In an effort to increase consistency, simplify maintenance and enable automated
change logs, we would like to kindly ask you to write _semantinc commit
messages_, as described in the [Conventional Commits
specification][conv-commits].

The general structure of _Conventional Commits_ is as follows:

```console
<type>[optional scope]: <description>

[optional body]

[optional footer]
```

Depending on the changes, please use one of the following **type** prefixes:

|    Type    |                                                                                                       Description                                                                                                             |
|:----------:|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|
| `build`    | The `build` type (formerly known as chore) is used to identify development changes related to the build system (involving scripts, configurations or tools) and package dependencies.                                         |
| `ci`       | The `ci` type is used to identify development changes related to the continuous integration and deployment system - involving scripts, configurations or tools.                                                               |
| `docs`     | The `docs` type is used to identify documentation changes related to the project - whether intended externally for the end users (in case of a library) or internally for the developers.                                     |
| `feat`     | The `feat` type is used to identify production changes related to new backward-compatible abilities or functionality.                                                                                                         |
| `fix`      | The `fix` type is used to identify production changes related to backward-compatible bug fixes.                                                                                                                               |
| `perf`     | The `perf` type is used to identify production changes related to backward-compatible performance improvements.                                                                                                               |
| `refactor` | The `refactor` type is used to identify development changes related to modifying the codebase, which neither adds a feature nor fixes a bug - such as removing redundant code, simplifying the code, renaming variables, etc. |
| `revert`   | For commits that revert one or more previous commits.                                                                                                                                                                         |
| `style`    | The `style`type is used to identify development changes related to styling the codebase, regardless of the meaning  such as indentations, semi-colons, quotes, trailing commas and so on.                                     |
| `test`     | The `test` type is used to identify development changes related to tests - such as refactoring existing tests or adding new tests.                                                                                            |

In order to ensure that the format of your commit message adheres to the
_Conventional Commits_ specification and the defined type vocabulary, you can
use the [dedicated linter][conv-commits-lint]. More information about
_Conventional Commits_ can also be found in this [blog
post][conv-commits-blog].

### Merging Your Code

Here is a check list that you can follow to make sure that code merges happen
smoothly:

1. [Open an issue](#submitting-issues) _first_ to give other contributors a
   chance to discuss the proposed changes (alternatively: assign yourself to
   one of the existing issues)
2. Clone the repository, create a feature branch off of the default branch
   (never commit changes to protected branches directly) and implement your
   code changes
3. If applicable, update relevant sections of the [main documentation][docs]
   and/or the [pipeline documentation][docs-pipeline]
4. Add or update tests; untested code will not be merged; refer to the
   [guidelines](#code-style-and-testing) above for details
5. Ensure that your coding style is in line with the
   [guidelines](#code-style-and-testing) described above
6. Ensure that all tests and linter checks configured in the
   [continuous integration][ci-cd] (CI) pipeline passes without issues
7. If necessary, clean up excessive commits with `git rebase`; cherry-pick and
   merge commits as you see fit; use concise and descriptive commit messages
8. Push your clean, tested and documented feature branch to the remote; make
   sure the [CI][ci-cd] pipeline passes
9. Issue a pull request against the default branch; follow the instructions in
   the [template][pull-request]; importantly, describe your changes in detail,
   yet with concise language, and do not forget to indicate which issue(s)
   the code changes resolve or refer to; assign a project maintainer to review
   your changes


## Becoming a Co-Maintainer

If you are as interested in the project as we are and have contributed some
code, suggested some features or reported bugs and have taken part in
discussions on where to go with the project, we will very likely want to have
you on board as a co-maintainer. If you are interested in that, please let us
know. You can reach us by [email][contact].


[bug-report]: .github/ISSUE_TEMPLATE/bug_report.md
[ci-cd]: .github/workflows/tests.yml
[coc]: CODE_OF_CONDUCT.md
[contact]: <mailto:zavolab-biozentrum@unibas.ch>
[conv-commits]: <https://www.conventionalcommits.org/en/v1.0.0-beta.2/#specification>
[conv-commits-blog]: <https://nitayneeman.com/posts/understanding-semantic-commit-messages-using-git-and-angular/>
[conv-commits-lint]: <https://github.com/conventional-changelog/commitlint>
[discussion-forum]: <https://github.com/zavolanlab/mirflowz/discussions>
[docs]: README.md
[docs-pipeline]: pipeline_documentation.md
[feat-request]: .github/ISSUE_TEMPLATE/feature_request.md
[git]: <https://git-scm.com/>
[git-flow]: <https://nvie.com/posts/a-successful-git-branching-model/>
[github]: <https://github.com>
[issue-tracker]: <https://github.com/zavolanlab/mirflowz/issues>
[pull-request]: PULL_REQUEST_TEMPLATE.md
[py]: <https://www.python.org/>
[py-flake8]: <https://flake8.pycqa.org>
[py-mypy]: <http://mypy-lang.org/>
[py-pylint]: <https://www.pylint.org/>
[py-pylint-conf]: pylint.cfg
[py-typing]: <https://docs.python.org/3/library/typing.html>
