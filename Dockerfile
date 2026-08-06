# Binder image for the definitive repository: FK-A, FK-B, and the shared model.
# Based on the official SageMath Binder environment.
FROM ghcr.io/sagemath/sage-binder-env:10.9

USER root

ARG NB_USER=user
ARG NB_UID=1000
ENV NB_USER=${NB_USER}
ENV NB_UID=${NB_UID}
ENV HOME=/home/${NB_USER}
ENV PATH="${HOME}/.local/bin:${PATH}"

# Where the base image actually puts Sage. It builds into /home/sage/sage and
# symlinks /usr/bin/sage, so `sage` is already on PATH; the tree is still needed
# by path for the bundled Jupyter. This is a build argument rather than an
# environment variable because SAGE_ROOT is meaningful to Sage itself, and this
# image should not be exporting one. Renaming uid 1000 below does not move the
# directory, so the path stays valid afterwards.
ARG SAGE_HOME=/home/sage/sage

# Binder assigns uid 1000. Reuse or rename the image's existing user so files
# mounted by repo2docker remain writable.
RUN if id -u ${NB_UID} >/dev/null 2>&1; then \
      EXISTING_USER=$(id -nu ${NB_UID}); \
      EXISTING_GROUP=$(id -gn ${NB_UID}); \
      if [ "${EXISTING_GROUP}" != "${NB_USER}" ]; then \
        groupmod -n ${NB_USER} ${EXISTING_GROUP} || true; \
      fi; \
      if [ "${EXISTING_USER}" != "${NB_USER}" ]; then \
        usermod -l ${NB_USER} -d ${HOME} ${EXISTING_USER}; \
        mkdir -p ${HOME}; \
        chown ${NB_USER}:${NB_USER} ${HOME}; \
      fi; \
    else \
      groupadd -g ${NB_UID} ${NB_USER} || true; \
      useradd -m -s /bin/bash -u ${NB_UID} -g ${NB_USER} ${NB_USER}; \
    fi

COPY --chown=${NB_USER}:${NB_USER} . ${HOME}/Hypergraph-Transversal-Thesis

USER ${NB_USER}
WORKDIR ${HOME}/Hypergraph-Transversal-Thesis

# No --user: Sage's Python is a virtualenv and rejects a user install. NB_USER
# is the image's own uid 1000, which owns the venv, so it can install into it.
RUN sage -python -m pip install --no-cache-dir -e ".[dev]" && \
    mkdir -p $(${SAGE_HOME}/venv/bin/jupyter --data-dir)/kernels && \
    ln -sf ${SAGE_HOME}/venv/share/jupyter/kernels/sagemath \
      $(${SAGE_HOME}/venv/bin/jupyter --data-dir)/kernels/sagemath
