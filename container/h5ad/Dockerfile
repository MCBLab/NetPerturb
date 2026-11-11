# Image for H5AD_TO_SEURAT, which turns an AnnData .h5ad (scanpy, CELLxGENE,
# ...) passed as --obj into the Seurat object the rest of the pipeline reads.
# modules/local/h5ad_to_seurat/main.nf names this tag, so build and push
# before the first non-stub run on an .h5ad:
#
#   docker build -t diegomscoelho/netperturb-h5ad:1.0 container/h5ad
#   docker push diegomscoelho/netperturb-h5ad:1.0
#
# Singularity pulls it from the same tag. If it cannot be pulled on your
# cluster, set the process's container to a local .sif in a -c config.
#
# Built on the very image DOWNSAMPLE and RANK_SCORE run in, pinned by digest,
# so the object is written by the same Seurat (4.4.0) and SeuratObject (5.0.0)
# that read it back: an object saved by a newer SeuratObject can carry slots
# or assay classes the older one cannot load. All that is added is hdf5r, to
# read the .h5ad directly -- no Python or anndata involved.
FROM juliaapolonio/scrank@sha256:d555e10af05a5e3469464f2e60ffe32114e87d9d0f7fb153ad334eee4bb2e86a

# hdf5r links against the system HDF5 library. Posit Package Manager has no
# Linux binary of it for the image's CRAN snapshot, so it is compiled from
# source, and its build applies a patch matching the HDF5 version it finds --
# hence `patch`, which the base image does not carry.
RUN apt-get update && apt-get install -y --no-install-recommends \
    libhdf5-dev \
    patch \
    && rm -rf /var/lib/apt/lists/*

# Same loud check as container/gsea and container/rquarto: fail the build
# rather than ship an image whose first task dies on a missing package.
RUN R -e "install.packages('hdf5r'); \
    if (!requireNamespace('hdf5r', quietly = TRUE)) stop('failed to install hdf5r')"

WORKDIR /workspace

CMD ["R"]
