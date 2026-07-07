FROM almalinux:9

# Install dependencies
RUN echo "#[INFO] Installing dependencies" && \
    dnf install -y python3 python3-pip git && \
    pip3 install mkdocs mkdocs-material pymdown-extensions plotly mkdocs-macros-plugin mkdocs-include-dir-to-nav mkdocs-include-markdown-plugin "markdown-exec[ansi]" && \
    dnf clean all && \
    pip3 cache purge


# Work directory
WORKDIR /app

# Copy files and launch
COPY . /app

# Command to launch the service
CMD ["mkdocs", "serve", "--livereload", "--dev-addr=0.0.0.0:8000"]