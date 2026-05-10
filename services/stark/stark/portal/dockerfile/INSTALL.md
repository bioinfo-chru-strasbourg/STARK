# Portal installation

## Python

```bash
conda create -n portail python=3.10
conda activate portail
pip install mkdocs mkdocs mkdocs-material pymdown-extensions plotly mkdocs-macros-plugin mkdocs-include-dir-to-nav mkdocs-include-markdown-plugin "markdown-exec[ansi]"
mkdocs serve --livereload
```
