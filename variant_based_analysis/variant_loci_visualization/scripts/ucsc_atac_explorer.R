suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
})

# ── Curated catalog of brain ATAC-seq track hubs ─────────────────────────────
#
# Each row = one bigwig track.  bw_url is either a full URL or a path
# relative to the hub base_url (resolved at fetch time).
#
# Sources:
#   BOCA  — Fullard et al. 2018, GSE96949 (NeuN+/NeuN- in 14 brain regions)
#   BOCA2 — Luo et al. 2020, Nat Commun  (Glut/GABA/Oligo/Mg-Astro in 3 regions)
#   AD    — Morabito et al. 2021, Nat Genet (pseudobulk snATAC in PFC, AD stages)
#
# Hint for the BOCA2 hub: visit labs.icahn.mssm.edu/roussos-lab/boca2/
# then click "UCSC Genome Browser" to get the exact hubUrl parameter.

ATAC_HUB_CATALOG <- data.frame(
  id = c(
    "boca_glut_ACC",   "boca_glut_dlPFC",  "boca_glut_V1",
    "boca_gaba_ACC",   "boca_gaba_dlPFC",  "boca_gaba_V1",
    "boca_oligo_ACC",  "boca_oligo_dlPFC", "boca_oligo_V1",
    "boca_mg_ACC",     "boca_mg_dlPFC",    "boca_mg_V1",
    "boca1_neuron_PFC","boca1_nonneuron_PFC",
    "ad_excN_PFC",     "ad_inhN_PFC",
    "ad_astro_PFC",    "ad_oligo_PFC",     "ad_micro_PFC",
    "encode_brain_cortex_atac"
  ),
  source = c(
    rep("BOCA2", 12),
    rep("BOCA1", 2),
    rep("AD_Regulome", 5),
    "ENCODE"
  ),
  cell_type = c(
    rep("Glutamatergic (Excitatory)", 3),
    rep("GABAergic (Inhibitory)", 3),
    rep("Oligodendrocyte", 3),
    rep("Microglia/Astrocyte", 3),
    "Neuron (NeuN+)", "Non-neuron (NeuN-)",
    "Excitatory neuron", "Inhibitory neuron",
    "Astrocyte", "Oligodendrocyte", "Microglia",
    "Brain cortex (bulk)"
  ),
  brain_region = c(
    "ACC", "dlPFC", "V1",
    "ACC", "dlPFC", "V1",
    "ACC", "dlPFC", "V1",
    "ACC", "dlPFC", "V1",
    "PFC", "PFC",
    "PFC", "PFC", "PFC", "PFC", "PFC",
    "Cortex"
  ),
  hub_url = c(
    # BOCA2 — parse hub.txt to get exact bigwig paths at runtime
    rep("https://bendlj01.u.hpc.mssm.edu/multireg2/hub.txt", 12),
    rep("https://bendlj01.u.hpc.mssm.edu/multireg/hub.txt", 2),
    rep("https://personal.broadinstitute.org/bjames/AD_snATAC/trackhub/hub.txt", 5),
    NA_character_   # ENCODE: resolved via API, not a hub
  ),
  # Known ENCODE ATAC-seq accession for brain cortex (hg38)
  # Source: ENCODE portal ENCSR526RCK (Homo sapiens brain cortex ATAC-seq)
  encode_accession = c(
    rep(NA_character_, 19),
    "ENCSR526RCK"
  ),
  # Short names used as track labels in Gviz
  label = c(
    "Glut ACC", "Glut dlPFC", "Glut V1",
    "GABA ACC", "GABA dlPFC", "GABA V1",
    "Oligo ACC", "Oligo dlPFC", "Oligo V1",
    "MgAstro ACC", "MgAstro dlPFC", "MgAstro V1",
    "Neuron PFC", "NonNeuron PFC",
    "ExcN PFC (AD)", "InhN PFC (AD)",
    "Astro PFC (AD)", "Oligo PFC (AD)", "Micro PFC (AD)",
    "Cortex (ENCODE)"
  ),
  recommended_for_excitatory = c(
    TRUE, TRUE, FALSE,    # ACC and dlPFC most relevant for GWAS loci
    FALSE, FALSE, FALSE,
    FALSE, FALSE, FALSE,
    FALSE, FALSE, FALSE,
    TRUE, FALSE,          # BOCA1 mixed neuron as fallback
    TRUE, FALSE, FALSE, FALSE, FALSE,
    TRUE                  # broad brain signal as sanity check
  ),
  stringsAsFactors = FALSE
)

# ── 1. UCSC REST API track browser ───────────────────────────────────────────

#' List tracks available in UCSC genome browser matching a pattern.
#'
#' @param pattern  Regex applied to track name AND shortLabel (case-insensitive)
#' @param genome   Assembly (default "hg38")
#' @return data.frame with columns: track, shortLabel, longLabel, type
browse_ucsc_tracks <- function(pattern = "atac", genome = "hg38") {
  url  <- sprintf("https://api.genome.ucsc.edu/list/tracks?genome=%s&trackLeavesOnly=1", genome)
  message("Querying UCSC API: ", url)
  resp <- tryCatch(
    GET(url, timeout(30)),
    error = function(e) { message("UCSC API unreachable: ", e$message); return(NULL) }
  )
  if (is.null(resp) || status_code(resp) != 200L) {
    message("UCSC API returned status ", status_code(resp)); return(NULL)
  }
  data <- fromJSON(rawToChar(resp$content), simplifyVector = FALSE)
  tracks_raw <- data[[genome]]
  if (is.null(tracks_raw)) { message("No tracks found for ", genome); return(NULL) }

  rows <- lapply(tracks_raw, function(t) {
    data.frame(
      track      = t$track       %||% NA_character_,
      shortLabel = t$shortLabel  %||% NA_character_,
      longLabel  = t$longLabel   %||% NA_character_,
      type       = t$type        %||% NA_character_,
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows)
  if (!is.null(pattern)) {
    hit <- grepl(pattern, df$track, ignore.case = TRUE) |
           grepl(pattern, df$shortLabel, ignore.case = TRUE) |
           grepl(pattern, df$longLabel, ignore.case = TRUE)
    df  <- df[hit, ]
  }
  if (nrow(df) == 0L) message("No native UCSC tracks matched pattern '", pattern, "'")
  df
}

# ── 2. Track hub parser ───────────────────────────────────────────────────────

#' Fetch and parse a UCSC track hub's trackDb.txt, returning a data.frame
#' of all bigwig signal tracks found.
#'
#' @param hub_url  URL to hub.txt (e.g. the hub_url column in ATAC_HUB_CATALOG)
#' @return data.frame: track, bigDataUrl, shortLabel, longLabel, type, abs_url
parse_hub_tracks <- function(hub_url) {
  # 1. Fetch hub.txt to find genomesFile
  message("Fetching hub: ", hub_url)
  hub_raw <- .fetch_text(hub_url)
  if (is.null(hub_raw)) return(NULL)

  genomes_file <- .hub_field(hub_raw, "genomesFile")
  if (is.null(genomes_file)) genomes_file <- "genomes.txt"
  base_url <- sub("/[^/]+$", "/", hub_url)

  # 2. Fetch genomes.txt to find hg38 trackDb
  genomes_url <- paste0(base_url, genomes_file)
  message("Fetching genomes file: ", genomes_url)
  genomes_raw <- .fetch_text(genomes_url)
  if (is.null(genomes_raw)) return(NULL)

  # Parse genomes file: find the hg38 block
  lines <- strsplit(genomes_raw, "\n")[[1L]]
  hg38_block_start <- which(trimws(lines) == "genome hg38")
  if (length(hg38_block_start) == 0L) {
    message("No hg38 block found in genomes file."); return(NULL)
  }
  # Extract trackDb line from hg38 block
  block_lines <- lines[hg38_block_start:min(hg38_block_start + 20L, length(lines))]
  trackdb_line <- grep("^trackDb", trimws(block_lines), value = TRUE)[1L]
  if (is.na(trackdb_line)) {
    message("No trackDb entry found for hg38"); return(NULL)
  }
  trackdb_file <- trimws(sub("^trackDb\\s+", "", trackdb_line))

  # 3. Fetch trackDb.txt — path is relative to hub base_url as written in genomes.txt
  trackdb_url <- paste0(base_url, trackdb_file)
  message("Fetching trackDb: ", trackdb_url)
  trackdb_raw <- .fetch_text(trackdb_url)
  if (is.null(trackdb_raw)) return(NULL)

  # 4. Parse trackDb stanzas; resolve relative bigDataUrls from the trackDb directory
  trackdb_dir <- dirname(trackdb_file)
  trackdb_base <- paste0(base_url, if (trackdb_dir == ".") "" else paste0(trackdb_dir, "/"))
  .parse_trackdb(trackdb_raw, base_url = trackdb_base)
}

# ── 3. Show catalog with filtering ───────────────────────────────────────────

#' Print a readable summary of available brain ATAC-seq tracks.
#'
#' @param filter_source  Keep only rows from this source (NULL = all)
#' @param excitatory_only  If TRUE, only show recommended_for_excitatory rows
#' @return data.frame (invisibly)
browse_catalog <- function(filter_source = NULL, excitatory_only = FALSE) {
  df <- ATAC_HUB_CATALOG
  if (!is.null(filter_source))
    df <- df[df$source %in% filter_source, ]
  if (excitatory_only)
    df <- df[df$recommended_for_excitatory, ]

  cat("\n── Brain ATAC-seq Track Catalog ─────────────────────────────────────────\n")
  cat(sprintf("%-25s %-8s %-30s %-18s %s\n",
              "ID", "Source", "Cell type", "Region", "Recommended(ExcN)"))
  cat(strrep("-", 95), "\n")
  for (i in seq_len(nrow(df))) {
    cat(sprintf("%-25s %-8s %-30s %-18s %s\n",
                df$id[i], df$source[i], df$cell_type[i],
                df$brain_region[i],
                if (df$recommended_for_excitatory[i]) "YES *" else ""))
  }
  cat("\nUse load_atac_hub_tracks(hub_url) to resolve bigwig URLs from a hub.\n")
  cat("Use build_catalog_track(id, chrom, from, to) to build a Gviz track.\n\n")
  invisible(df)
}

# ── 4. Gviz track builder from bigwig URL ─────────────────────────────────────

#' Build a Gviz DataTrack from a bigwig URL (local file or http/https).
#' Requires rtracklayer.  Returns NULL silently on failure.
#'
#' @param bw_url     Path or URL to a bigwig file
#' @param chrom      UCSC chromosome, e.g. "chr5"
#' @param from       Start coordinate (1-based, as used by Gviz)
#' @param to         End coordinate (1-based)
#' @param name       Track label (default derived from URL)
#' @param color      Bar color (default brain blue)
#' @param track_type Gviz type arg passed to DataTrack (default "histogram")
#' @return DataTrack or NULL
build_bw_track <- function(bw_url, chrom, from, to,
                            name       = NULL,
                            color      = "#1B6CA8",
                            track_type = "histogram",
                            genome     = "hg38") {
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    message("rtracklayer required: BiocManager::install('rtracklayer')"); return(NULL)
  }
  if (is.null(name))
    name <- sub("\\.[^.]+$", "", basename(bw_url))

  which_gr <- GenomicRanges::GRanges(chrom, IRanges::IRanges(from, to))

  message("  Importing bigwig: ", name, "  [", chrom, ":", from, "-", to, "]")
  bw_data <- tryCatch(
    rtracklayer::import.bw(bw_url, which = which_gr),
    error = function(e) { message("  [WARN] BigWig import failed: ", e$message); NULL }
  )
  if (is.null(bw_data) || length(bw_data) == 0L) {
    message("  [WARN] No data in region for: ", name); return(NULL)
  }

  Gviz::DataTrack(
    range            = bw_data,
    type             = track_type,
    genome           = genome,
    name             = name,
    col.histogram    = color,
    fill.histogram   = color,
    col              = color,
    background.title = "white",
    col.title        = "black"
  )
}

#' Build a Gviz track from the catalog by track id.
#' Resolves the bigwig URL by parsing the hub's trackDb.txt at runtime.
#'
#' @param catalog_id  One of ATAC_HUB_CATALOG$id
#' @param chrom, from, to  Locus as used by Gviz (1-based coords)
#' @return DataTrack or NULL
build_catalog_track <- function(catalog_id, chrom, from, to, ...) {
  row <- ATAC_HUB_CATALOG[ATAC_HUB_CATALOG$id == catalog_id, ]
  if (nrow(row) == 0L) { message("Unknown catalog id: ", catalog_id); return(NULL) }

  hub_url <- row$hub_url
  if (is.na(hub_url)) {
    message("No hub_url for ", catalog_id, " — resolve manually."); return(NULL)
  }

  tracks_df <- parse_hub_tracks(hub_url)
  if (is.null(tracks_df)) return(NULL)

  # Match by cell type keywords in shortLabel / longLabel
  kw <- .id_to_keywords(catalog_id)
  idx <- Reduce(`&`, lapply(kw, function(k)
    grepl(k, tracks_df$shortLabel, ignore.case = TRUE) |
    grepl(k, tracks_df$longLabel,  ignore.case = TRUE)))

  bw_rows <- tracks_df[idx & grepl("bigwig|bigWig|bw", tracks_df$type, ignore.case = TRUE), ]
  if (nrow(bw_rows) == 0L) {
    message("No matching bigwig track for id '", catalog_id, "' in hub. ",
            "Inspect parse_hub_tracks('", hub_url, "') manually.")
    return(NULL)
  }
  if (nrow(bw_rows) > 1L) {
    message("Multiple matches for '", catalog_id, "' — using first: ", bw_rows$shortLabel[1L])
    bw_rows <- bw_rows[1L, ]
  }
  build_bw_track(bw_rows$abs_url, chrom, from, to,
                 name  = row$label,
                 color = .source_color(row$source), ...)
}

# ── 5. Load and display hub tracks (sanity-check helper) ──────────────────────

#' Fetch and list all bigwig tracks from a hub for hg38.
#' Use this to browse what's actually in a hub before picking a track ID.
#'
#' @param hub_url  URL to hub.txt
#' @return data.frame with track metadata + resolved bigwig URL
load_atac_hub_tracks <- function(hub_url) {
  df <- parse_hub_tracks(hub_url)
  if (is.null(df)) return(invisible(NULL))
  bw_df <- df[grepl("bigwig|bw", df$type, ignore.case = TRUE) & !is.na(df$bigDataUrl), ]
  if (nrow(bw_df) == 0L) { message("No bigwig tracks found in hub."); return(invisible(df)) }
  cat("\n── Bigwig tracks in hub ─────────────────────────────────────────────────\n")
  cat(sprintf("%-35s %-35s %s\n", "shortLabel", "longLabel", "abs_url"))
  cat(strrep("-", 100), "\n")
  for (i in seq_len(nrow(bw_df))) {
    cat(sprintf("%-35s %-35s %s\n",
                substr(bw_df$shortLabel[i], 1, 33),
                substr(bw_df$longLabel[i],  1, 33),
                bw_df$abs_url[i]))
  }
  invisible(bw_df)
}

# ── 6. UcscTrack wrapper (for UCSC-native bigwig tracks) ──────────────────────

#' Build a Gviz UcscTrack for a track that exists natively on genome.ucsc.edu.
#' This only works for tracks in the UCSC MySQL tables (not external hubs).
#' For track hubs, use build_bw_track() instead.
#'
#' @param track  UCSC track table name, e.g. "encRegTfbsClusters"
#' @param genome, chrom, from, to  Locus
#' @param type   Gviz track type (default "bw")
#' @param name   Track label
build_ucsc_native_track <- function(track, genome = "hg38", chrom, from, to,
                                     type = "bw", name = track, ...) {
  tryCatch(
    Gviz::UcscTrack(
      genome           = genome,
      chromosome       = chrom,
      track            = track,
      from             = from,
      to               = to,
      trackType        = type,
      name             = name,
      background.title = "white",
      col.title        = "black",
      ...
    ),
    error = function(e) {
      message("UcscTrack failed for '", track, "': ", e$message)
      NULL
    }
  )
}

# ── 7. ENCODE REST API — resolve bigwig URL from experiment accession ─────────

#' Resolve the download URL of a released hg38 bigWig from an ENCODE experiment.
#'
#' Queries the ENCODE portal JSON API and returns the first released bigWig
#' file whose output_type matches \code{output_type_pattern}.
#'
#' @param accession  ENCODE experiment accession, e.g. "ENCSR526RCK"
#' @param output_type_pattern  Regex matched against output_type field
#'                             (default matches "signal p-value" or "signal")
#' @param assembly   Genome assembly string as used by ENCODE (default "GRCh38")
#' @return Character URL or NULL on failure
resolve_encode_bw_url <- function(accession,
                                   output_type_pattern = "signal",
                                   assembly = "GRCh38") {
  url  <- sprintf("https://www.encodeproject.org/experiments/%s/?format=json",
                  accession)
  message("Querying ENCODE portal: ", url)
  resp <- tryCatch(
    GET(url, timeout(30)),
    error = function(e) { message("ENCODE API unreachable: ", e$message); NULL }
  )
  if (is.null(resp) || status_code(resp) != 200L) {
    message("ENCODE API status: ", if (!is.null(resp)) status_code(resp) else "N/A")
    return(NULL)
  }
  data  <- fromJSON(rawToChar(resp$content), simplifyVector = FALSE)
  files <- data$files
  if (is.null(files) || length(files) == 0L) {
    message("No files listed for ENCODE experiment ", accession); return(NULL)
  }

  # Prefer released bigWig for target assembly matching output_type_pattern
  .keep <- function(f, require_otype = TRUE) {
    isTRUE(f[["file_type"]] == "bigWig") &&
    isTRUE(grepl(assembly, f[["assembly"]] %||% "", ignore.case = TRUE)) &&
    isTRUE((f[["status"]] %||% "") == "released") &&
    (!require_otype || grepl(output_type_pattern, f[["output_type"]] %||% "",
                             ignore.case = TRUE))
  }
  bw <- Filter(function(f) .keep(f, TRUE),  files)
  if (length(bw) == 0L)
    bw <- Filter(function(f) .keep(f, FALSE), files)   # relax output_type filter
  if (length(bw) == 0L) {
    message("No released GRCh38 bigWig found for ", accession); return(NULL)
  }
  href <- bw[[1L]][["href"]]
  if (is.null(href)) { message("No href for bigWig file"); return(NULL) }
  paste0("https://www.encodeproject.org", href)
}

# ── Internals ─────────────────────────────────────────────────────────────────

`%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && nchar(a) > 0) a else b

.fetch_text <- function(url, timeout_s = 20L) {
  tryCatch({
    resp <- GET(url, timeout(timeout_s))
    if (status_code(resp) != 200L) { message("HTTP ", status_code(resp), " for ", url); return(NULL) }
    rawToChar(resp$content)
  }, error = function(e) { message("Fetch error: ", e$message); NULL })
}

.hub_field <- function(text, field) {
  m <- regmatches(text, regexpr(paste0("(?m)^", field, "\\s+(\\S+)"), text, perl = TRUE))
  if (length(m) == 0L) return(NULL)
  trimws(sub(paste0("^", field, "\\s+"), "", m))
}

.parse_trackdb <- function(text, base_url = "") {
  # Split into stanzas (blank-line-separated)
  stanzas <- strsplit(text, "\n\\s*\n")[[1L]]
  rows <- lapply(stanzas, function(s) {
    track      <- .stanza_field(s, "track")
    bdu        <- .stanza_field(s, "bigDataUrl")
    short      <- .stanza_field(s, "shortLabel")
    long       <- .stanza_field(s, "longLabel")
    type       <- .stanza_field(s, "type")
    if (is.null(track)) return(NULL)
    # Resolve relative URL
    abs_url <- if (!is.null(bdu) && grepl("^https?://", bdu)) {
      bdu
    } else if (!is.null(bdu)) {
      paste0(base_url, bdu)
    } else NA_character_
    data.frame(track = track, bigDataUrl = bdu %||% NA_character_,
               shortLabel = short %||% NA_character_,
               longLabel  = long  %||% NA_character_,
               type       = type  %||% NA_character_,
               abs_url    = abs_url,
               stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) return(NULL)
  do.call(rbind, rows)
}

.stanza_field <- function(stanza, field) {
  m <- regmatches(stanza, regexpr(paste0("(?m)^\\s*", field, "\\s+(.+)"), stanza, perl = TRUE))
  if (length(m) == 0L) return(NULL)
  trimws(sub(paste0("^\\s*", field, "\\s+"), "", m))
}

# Map catalog id → keywords for matching in trackDb shortLabel/longLabel
.id_to_keywords <- function(id) {
  map <- list(
    boca_glut_ACC    = c("glutamatergic", "ACC"),
    boca_glut_dlPFC  = c("glutamatergic", "dlPFC"),
    boca_glut_V1     = c("glutamatergic", "V1"),
    boca_gaba_ACC    = c("GABAergic", "ACC"),
    boca_gaba_dlPFC  = c("GABAergic", "dlPFC"),
    boca_gaba_V1     = c("GABAergic", "V1"),
    boca_oligo_ACC   = c("oligodendrocyte", "ACC"),
    boca_oligo_dlPFC = c("oligodendrocyte", "dlPFC"),
    boca_oligo_V1    = c("oligodendrocyte", "V1"),
    boca_mg_ACC      = c("microglia", "ACC"),
    boca_mg_dlPFC    = c("microglia", "dlPFC"),
    boca_mg_V1       = c("microglia", "V1"),
    boca1_neuron_PFC    = c("NeuN"),
    boca1_nonneuron_PFC = c("NeuN-"),
    ad_excN_PFC      = c("excitatory"),
    ad_inhN_PFC      = c("inhibitory"),
    ad_astro_PFC     = c("astrocyte"),
    ad_oligo_PFC     = c("oligodendrocyte"),
    ad_micro_PFC     = c("microglia")
  )
  map[[id]] %||% list(id)
}

.source_color <- function(source) {
  switch(source,
    BOCA2       = "#1B6CA8",
    BOCA1       = "#7EB8D9",
    AD_Regulome = "#D45E1A",
    ENCODE      = "#444444",
    "#666666"
  )
}

# ── 8. UCSC URL parser and convenience entry points ──────────────────────────

#' Parse query parameters from a URL string into a named list.
.parse_url_params <- function(url) {
  query_str <- sub("^[^?]*\\?", "", url)
  pairs     <- strsplit(query_str, "&")[[1L]]
  names_    <- sub("=.*$", "", pairs)
  vals      <- vapply(pairs, function(p) utils::URLdecode(sub("^[^=]+=", "", p)),
                      character(1L))
  setNames(as.list(vals), names_)
}

#' Parse a UCSC Genome Browser URL into genome, hub_url, and region components.
#'
#' @param ucsc_url  Full browser URL, e.g. from the address bar or a hub link
#' @return Named list with \code{genome}, \code{hub_url}, \code{chrom}, \code{from}, \code{to}
parse_ucsc_url <- function(ucsc_url) {
  p      <- .parse_url_params(ucsc_url)
  genome <- p[["db"]] %||% "hg38"

  hub_url <- p[["hubUrl"]]

  region <- NULL
  pos    <- p[["position"]]
  if (!is.null(pos) && nchar(pos) > 0L) {
    m <- regmatches(pos, regexec("^([^:]+):(\\d+)-(\\d+)$", pos))[[1L]]
    if (length(m) == 4L)
      region <- list(chrom = m[2L], from = as.integer(m[3L]), to = as.integer(m[4L]))
  }

  list(genome = genome, hub_url = hub_url,
       chrom  = region$chrom, from = region$from, to = region$to)
}

#' Build Gviz DataTracks for every bigwig in a track hub over a genomic region.
#'
#' @param hub_url    URL to hub.txt
#' @param chrom, from, to  Locus (1-based, UCSC-style)
#' @param genome     Assembly string passed to Gviz (default "hg38")
#' @param track_type Gviz type argument (default "histogram")
#' @param colors     Optional character vector of colors (recycled)
#' @return List of DataTrack objects (NULLs dropped)
build_hub_tracks <- function(hub_url, chrom, from, to,
                              genome     = "hg38",
                              track_type = "histogram",
                              colors     = c("#1B6CA8","#D45E1A","#2E8B57",
                                             "#8B2E8B","#8B7D2E","#666666")) {
  tracks_df <- parse_hub_tracks(hub_url)
  if (is.null(tracks_df)) return(list())

  bw_df <- tracks_df[grepl("bigwig|bw", tracks_df$type, ignore.case = TRUE) &
                     !is.na(tracks_df$bigDataUrl), ]
  if (nrow(bw_df) == 0L) { message("No bigwig tracks found in hub."); return(list()) }

  message("Building ", nrow(bw_df), " bigwig track(s) for ",
          chrom, ":", from, "-", to)

  track_list <- lapply(seq_len(nrow(bw_df)), function(i) {
    col <- colors[((i - 1L) %% length(colors)) + 1L]
    build_bw_track(bw_df$abs_url[i], chrom, from, to,
                   name       = bw_df$shortLabel[i],
                   color      = col,
                   track_type = track_type,
                   genome     = genome)
  })
  Filter(Negate(is.null), track_list)
}

#' One-shot: parse a UCSC Genome Browser URL and return ready-to-use Gviz tracks.
#'
#' Example:
#'   tracks <- build_tracks_from_ucsc_url(
#'     "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&hubUrl=https://ggoma.s3.amazonaws.com/hub_basic.txt&position=chr19:35900492-35912218"
#'   )
#'   library(Gviz)
#'   plotTracks(c(GenomeAxisTrack(), tracks))
#'
#' @param ucsc_url  Full UCSC Genome Browser URL containing hubUrl= and position=
#' @param ...       Extra arguments forwarded to \code{build_hub_tracks}
#' @return List of DataTrack objects
build_tracks_from_ucsc_url <- function(ucsc_url, ...) {
  parsed <- parse_ucsc_url(ucsc_url)
  if (is.null(parsed$hub_url)) { message("No hubUrl= parameter found in URL"); return(list()) }
  if (is.null(parsed$chrom))   { message("No position= parameter found in URL"); return(list()) }
  message("Hub    : ", parsed$hub_url)
  message("Region : ", parsed$chrom, ":", parsed$from, "-", parsed$to)
  build_hub_tracks(parsed$hub_url, parsed$chrom, parsed$from, parsed$to,
                   genome = parsed$genome, ...)
}

# ── Startup message ───────────────────────────────────────────────────────────
message("\n── ucsc_atac_explorer.R loaded ─────────────────────────────────────────────")
message("Quick-start:")
message("  # Fastest path — paste any UCSC browser URL directly:")
message("  tracks <- build_tracks_from_ucsc_url('http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&hubUrl=https://ggoma.s3.amazonaws.com/hub_basic.txt&position=chr19:35900492-35912218')")
message("  plotTracks(c(GenomeAxisTrack(), tracks))")
message("")
message("  # Or build from explicit hub URL + region:")
message("  tracks <- build_hub_tracks('https://ggoma.s3.amazonaws.com/hub_basic.txt', 'chr19', 35900492, 35912218)")
message("")
message("  browse_catalog()                         # show curated catalog")
message("  browse_catalog(excitatory_only = TRUE)   # excitatory-neuron tracks only")
message("  load_atac_hub_tracks('https://ggoma.s3.amazonaws.com/hub_basic.txt')  # list bigwigs in hub")
message("  build_catalog_track('boca_glut_dlPFC', chrom, from, to)  # catalog → Gviz track")
message("  resolve_encode_bw_url('ENCSR526RCK')                     # ENCODE API → bigwig URL")
