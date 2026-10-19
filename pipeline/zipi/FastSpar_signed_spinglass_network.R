# ============================================================
# 有符号微生物网络分析脚本：导入 FastSpar 相关性与 P 值矩阵
#
# 核心变化：
# 1. 保留相关系数原始正负号：E(g)$correlation
# 2. 保存边类型：E(g)$sign / E(g)$label
# 3. E(g)$weight 仅保存 abs(correlation)，供绘图和不支持负权重的函数使用
# 4. 使用 cluster_spinglass(..., implementation = "neg") 进行有符号模块划分
# 5. 对每个连通分量分别运行 signed spin-glass，再统一重编号为 M1、M2...
# 6. Zi/Pi 仍按 microeco 的标准非加权定义计算，但模块来源改为有符号模块
# 7. 节点表额外输出正/负、模块内/模块间连接数
# 8. 仅对 network_plot.html 执行自适应展示稀疏化：
#    平均度 > 12 时，每个节点按 |r| 保留最强 3 条边并取并集
#    完整网络仍用于模块、Zi/Pi、网络指标、文本结果和 GEXF
# 9. network_plot.html 中正相关边与负相关边使用不同颜色，并显示边符号图例
# 10. 同时输出正相关、负相关以及正负纵向双面板 coabundance HTML
# 11. 输出正负 Phylum 汇总矩阵、配对汇总和差异诊断，确认两图是否真正不同
# 12. 增大 chord 图边距并清理 p__ 前缀，避免标签被裁切
# ============================================================

library(tidyverse)
library(microeco)
library(igraph)
library(htmlwidgets)
library(networkD3)
library(gridExtra)
library(optparse)

message(
    "运行环境：R ", getRversion(),
    "；igraph ", as.character(utils::packageVersion("igraph")),
    "；microeco ", as.character(utils::packageVersion("microeco"))
)

check_signed_spinglass_support <- function() {
    test_mat <- matrix(
        c(
             0,  0.8, -0.8,
             0.8, 0,  -0.7,
            -0.8, -0.7, 0
        ),
        nrow = 3,
        byrow = TRUE
    )
    test_g <- igraph::graph_from_adjacency_matrix(
        test_mat,
        mode = "undirected",
        weighted = TRUE,
        diag = FALSE
    )

    test_fit <- tryCatch(
        igraph::cluster_spinglass(
            test_g,
            weights = E(test_g)$weight,
            spins = 3,
            parupdate = FALSE,
            update.rule = "config",
            gamma = 1,
            implementation = "neg",
            gamma.minus = 1
        ),
        error = identity
    )

    if (inherits(test_fit, "error")) {
        err_print <- paste(capture.output(print(test_fit)), collapse = " | ")
        err_str <- paste(capture.output(str(test_fit)), collapse = " | ")
        stop(
            paste0(
                "当前 igraph 无法通过最小 signed spin-glass 自检。",
                " igraph版本=", as.character(utils::packageVersion("igraph")),
                "；conditionMessage=", conditionMessage(test_fit),
                "；print(error)=", err_print,
                "；str(error)=", err_str,
                "。请先升级/重装 igraph 后再运行。"
            ),
            call. = FALSE
        )
    }

    message(
        "signed spin-glass 最小自检通过；测试模块数=",
        length(unique(igraph::membership(test_fit)))
    )
}
check_signed_spinglass_support()

empty_plot <- function(label) {
    ggplot() +
        annotate("text", x = 0, y = 0, label = label, size = 6) +
        theme_void()
}

write_empty_html <- function(file, label) {
    writeLines(
        c(
            "<!doctype html>",
            "<html><head><meta charset=\"utf-8\"><title>No plot</title></head>",
            paste0("<body><p>", label, "</p></body></html>")
        ),
        con = file
    )
}

safe_normalize_dir <- function(path) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    normalizePath(path, mustWork = TRUE)
}

build_basic_node_table <- function(g, dataset) {
    node_names <- V(g)$name

    safe_metric <- function(x) {
        out <- tryCatch(x, error = function(e) numeric())
        if (length(out) != length(node_names)) {
            rep(NA_real_, length(node_names))
        } else {
            as.numeric(out)
        }
    }

    degree_vals <- safe_metric(igraph::degree(g))
    betweenness_vals <- safe_metric(igraph::betweenness(g, directed = FALSE))
    closeness_vals <- safe_metric(
        igraph::closeness(g, mode = "all", normalized = FALSE, weights = NA)
    )
    eig <- safe_metric(
        igraph::eigen_centrality(
            g,
            directed = FALSE,
            weights = E(g)$weight
        )$vector
    )
    abundance_vals <- safe_metric(V(g)$Abundance)

    module_vals <- as.character(V(g)$module)
    if (length(module_vals) != length(node_names)) {
        module_vals <- rep(NA_character_, length(node_names))
    }

    out <- data.frame(
        name = node_names,
        degree = degree_vals,
        betweenness_centrality = betweenness_vals,
        closeness_centrality = closeness_vals,
        eigenvector_centrality = eig,
        Abundance = abundance_vals,
        module = module_vals,
        z = NA_real_,
        p = NA_real_,
        taxa_roles = NA_character_,
        stringsAsFactors = FALSE
    )

    tax_table <- dataset$tax_table
    tax_match <- tax_table[node_names, , drop = FALSE]
    tax_match <- tax_match[, setdiff(colnames(tax_match), colnames(out)), drop = FALSE]
    cbind(out, tax_match)
}

# ------------------------------------------------------------------
# 对有符号网络逐个连通分量运行 signed spin-glass
# ------------------------------------------------------------------
assign_signed_spinglass_modules <- function(
    g,
    spins = 25L,
    gamma = 1,
    gamma_minus = 1,
    seed = 12345L
) {
    if (igraph::is_directed(g)) {
        stop("signed spin-glass 模块划分要求本脚本输入无向网络。")
    }
    if (!"correlation" %in% igraph::edge_attr_names(g)) {
        stop("网络缺少 E(g)$correlation，无法进行有符号模块划分。")
    }
    if (any(!is.finite(E(g)$correlation))) {
        stop("E(g)$correlation 中存在非有限值。")
    }

    comp <- igraph::components(g, mode = "weak")
    component_ids <- sort(unique(comp$membership))

    raw_module <- setNames(rep(NA_integer_, igraph::vcount(g)), V(g)$name)
    component_records <- vector("list", length(component_ids))
    next_global_module <- 1L

    for (ii in seq_along(component_ids)) {
        component_id <- component_ids[ii]
        vertex_index <- which(comp$membership == component_id)
        subg <- igraph::induced_subgraph(g, vids = vertex_index)

        n_nodes <- igraph::vcount(subg)
        n_edges <- igraph::ecount(subg)
        n_positive <- sum(E(subg)$correlation > 0)
        n_negative <- sum(E(subg)$correlation < 0)
        component_seed <- as.integer(seed) + as.integer(component_id) - 1L

        message(
            "signed spin-glass：连通分量 ", component_id,
            "，节点=", n_nodes,
            "，边=", n_edges,
            "，正边=", n_positive,
            "，负边=", n_negative
        )

        # 删除孤立节点后原则上不会出现单节点分量，但仍保留安全处理。
        if (n_nodes == 1L || n_edges == 0L) {
            local_membership <- setNames(1L, V(subg)$name)
            component_method <- "singleton_fallback"
        } else if (n_nodes == 2L) {
            # 两节点正相关：同一模块；两节点负相关：分别作为两个模块。
            if (E(subg)$correlation[1] > 0) {
                local_membership <- setNames(c(1L, 1L), V(subg)$name)
                component_method <- "two_node_positive_same_module"
            } else {
                local_membership <- setNames(c(1L, 2L), V(subg)$name)
                component_method <- "two_node_negative_separate_modules"
            }
        } else {
            local_spins <- max(2L, min(as.integer(spins), n_nodes))
            set.seed(component_seed)

            message(
                "分量检查：connected=", igraph::is_connected(subg),
                "；simple=", igraph::is_simple(subg),
                "；weight范围=[",
                signif(min(E(subg)$correlation), 6), ", ",
                signif(max(E(subg)$correlation), 6), "]",
                "；零权重数=", sum(E(subg)$correlation == 0),
                "；非有限权重数=", sum(!is.finite(E(subg)$correlation))
            )

            fit <- tryCatch(
                igraph::cluster_spinglass(
                    subg,
                    weights = E(subg)$correlation,
                    spins = local_spins,
                    parupdate = FALSE,
                    update.rule = "config",
                    gamma = as.numeric(gamma),
                    implementation = "neg",
                    gamma.minus = as.numeric(gamma_minus)
                ),
                error = identity
            )

            if (inherits(fit, "error")) {
                err_print <- paste(capture.output(print(fit)), collapse = " | ")
                err_str <- paste(capture.output(str(fit)), collapse = " | ")
                stop(
                    paste0(
                        "连通分量 ", component_id,
                        " 的 signed spin-glass 划分失败。",
                        " igraph版本=", as.character(utils::packageVersion("igraph")),
                        "；conditionMessage=", conditionMessage(fit),
                        "；geterrmessage=", geterrmessage(),
                        "；print(error)=", err_print,
                        "；str(error)=", err_str
                    ),
                    call. = FALSE
                )
            }

            local_membership <- igraph::membership(fit)
            local_membership <- setNames(
                as.integer(local_membership),
                V(subg)$name
            )
            component_method <- "signed_spinglass_neg"
        }

        local_levels <- sort(unique(local_membership))
        local_to_global <- setNames(
            seq.int(
                from = next_global_module,
                length.out = length(local_levels)
            ),
            as.character(local_levels)
        )

        raw_module[names(local_membership)] <- unname(
            local_to_global[as.character(local_membership)]
        )
        next_global_module <- next_global_module + length(local_levels)

        component_records[[ii]] <- data.frame(
            component = component_id,
            nodes = n_nodes,
            edges = n_edges,
            positive_edges = n_positive,
            negative_edges = n_negative,
            modules = length(local_levels),
            method = component_method,
            seed = component_seed,
            spins = if (n_nodes > 2L) max(2L, min(as.integer(spins), n_nodes)) else NA_integer_,
            gamma = as.numeric(gamma),
            gamma_minus = as.numeric(gamma_minus),
            stringsAsFactors = FALSE
        )
    }

    if (any(is.na(raw_module))) {
        stop("部分节点未获得模块编号，停止分析。")
    }

    # 按模块节点数从大到小统一命名 M1、M2...；同大小时按最小节点名排序。
    raw_module_character <- as.character(raw_module)
    module_size <- table(raw_module_character)
    module_first_name <- tapply(
        names(raw_module),
        raw_module_character,
        function(x) sort(x)[1]
    )

    module_order <- names(module_size)[
        order(
            -as.integer(module_size),
            module_first_name[names(module_size)]
        )
    ]
    module_map <- setNames(
        paste0("M", seq_along(module_order)),
        module_order
    )

    final_module <- unname(module_map[raw_module_character])
    names(final_module) <- names(raw_module)

    g <- igraph::set_vertex_attr(
        g,
        "module",
        value = final_module[V(g)$name]
    )

    g <- igraph::set_graph_attr(g, "module_method", "signed_spinglass_neg_componentwise")
    g <- igraph::set_graph_attr(g, "module_seed", as.integer(seed))
    g <- igraph::set_graph_attr(g, "module_spins", as.integer(spins))
    g <- igraph::set_graph_attr(g, "module_gamma", as.numeric(gamma))
    g <- igraph::set_graph_attr(g, "module_gamma_minus", as.numeric(gamma_minus))

    list(
        graph = g,
        component_summary = dplyr::bind_rows(component_records)
    )
}

# ------------------------------------------------------------------
# 给节点表增加正负边和模块内/模块间边的计数
# ------------------------------------------------------------------
calculate_signed_node_metrics <- function(g) {
    node_names <- V(g)$name
    module_map <- setNames(as.character(V(g)$module), node_names)

    out <- lapply(node_names, function(node_name) {
        edge_ids <- igraph::incident(g, node_name, mode = "all")

        if (length(edge_ids) == 0L) {
            return(data.frame(
                name = node_name,
                positive_degree = 0L,
                negative_degree = 0L,
                within_positive_degree = 0L,
                within_negative_degree = 0L,
                between_positive_degree = 0L,
                between_negative_degree = 0L,
                positive_strength = 0,
                negative_strength_abs = 0,
                stringsAsFactors = FALSE
            ))
        }

        correlations <- E(g)[edge_ids]$correlation
        edge_ends <- igraph::ends(g, edge_ids, names = TRUE)
        other_nodes <- ifelse(
            edge_ends[, 1] == node_name,
            edge_ends[, 2],
            edge_ends[, 1]
        )
        same_module <- module_map[other_nodes] == module_map[node_name]

        data.frame(
            name = node_name,
            positive_degree = sum(correlations > 0),
            negative_degree = sum(correlations < 0),
            within_positive_degree = sum(correlations > 0 & same_module),
            within_negative_degree = sum(correlations < 0 & same_module),
            between_positive_degree = sum(correlations > 0 & !same_module),
            between_negative_degree = sum(correlations < 0 & !same_module),
            positive_strength = sum(correlations[correlations > 0]),
            negative_strength_abs = sum(abs(correlations[correlations < 0])),
            stringsAsFactors = FALSE
        )
    })

    dplyr::bind_rows(out)
}

# ------------------------------------------------------------------
# 输出边表，并标明是否为模块内部边
# ------------------------------------------------------------------
build_signed_edge_table <- function(g) {
    edge_table <- igraph::as_data_frame(g, what = "edges")
    module_map <- setNames(as.character(V(g)$module), V(g)$name)

    edge_table$from_module <- unname(module_map[edge_table$from])
    edge_table$to_module <- unname(module_map[edge_table$to])
    edge_table$same_module <- edge_table$from_module == edge_table$to_module
    edge_table
}

# ------------------------------------------------------------------
# 按模块汇总内部及外部正负边，便于测试有符号模块是否合理
# ------------------------------------------------------------------
build_module_summary <- function(g, edge_table = NULL) {
    if (is.null(edge_table)) {
        edge_table <- build_signed_edge_table(g)
    }

    module_names <- sort(
        unique(as.character(V(g)$module)),
        method = "radix"
    )

    out <- lapply(module_names, function(module_name) {
        module_nodes <- V(g)$name[V(g)$module == module_name]

        internal <- edge_table[
            edge_table$from_module == module_name &
                edge_table$to_module == module_name,
            ,
            drop = FALSE
        ]

        external <- edge_table[
            xor(
                edge_table$from_module == module_name,
                edge_table$to_module == module_name
            ),
            ,
            drop = FALSE
        ]

        internal_total <- nrow(internal)
        internal_positive <- sum(internal$correlation > 0)
        internal_negative <- sum(internal$correlation < 0)
        external_total <- nrow(external)
        external_positive <- sum(external$correlation > 0)
        external_negative <- sum(external$correlation < 0)

        data.frame(
            module = module_name,
            nodes = length(module_nodes),
            internal_edges = internal_total,
            internal_positive_edges = internal_positive,
            internal_negative_edges = internal_negative,
            internal_positive_fraction = if (internal_total > 0) internal_positive / internal_total else NA_real_,
            external_edges = external_total,
            external_positive_edges = external_positive,
            external_negative_edges = external_negative,
            external_negative_fraction = if (external_total > 0) external_negative / external_total else NA_real_,
            stringsAsFactors = FALSE
        )
    })

    summary_table <- dplyr::bind_rows(out)
    module_number <- suppressWarnings(as.integer(sub("^M", "", summary_table$module)))
    summary_table[order(module_number), , drop = FALSE]
}


# ------------------------------------------------------------------
# 仅用于网页可视化的自适应边稀疏化
#
# 规则：
# 1. 完整网络平均度 <= avg_degree_threshold：显示全部边；
# 2. 完整网络平均度 > avg_degree_threshold：
#    对每个节点按 |correlation| 从大到小保留 top_k 条边，
#    再对所有节点选中的边取并集并去重；
# 3. 不删除任何节点；仅返回绘图网络副本，不改变完整分析网络。
# ------------------------------------------------------------------
build_plot_network <- function(
    g,
    avg_degree_threshold = 12,
    top_k = 3L
) {
    if (igraph::vcount(g) == 0L) {
        return(list(
            graph = g,
            applied = FALSE,
            original_nodes = 0L,
            original_edges = 0L,
            original_average_degree = NA_real_,
            plot_nodes = 0L,
            plot_edges = 0L,
            plot_average_degree = NA_real_,
            avg_degree_threshold = avg_degree_threshold,
            top_k = top_k
        ))
    }

    original_average_degree <- mean(igraph::degree(g))

    if (
        !is.finite(original_average_degree) ||
        original_average_degree <= avg_degree_threshold ||
        igraph::ecount(g) == 0L
    ) {
        return(list(
            graph = g,
            applied = FALSE,
            original_nodes = igraph::vcount(g),
            original_edges = igraph::ecount(g),
            original_average_degree = original_average_degree,
            plot_nodes = igraph::vcount(g),
            plot_edges = igraph::ecount(g),
            plot_average_degree = original_average_degree,
            avg_degree_threshold = avg_degree_threshold,
            top_k = top_k
        ))
    }

    if (top_k < 1L) {
        stop("用于网页绘图的 top_k 必须不小于 1。")
    }

    if ("correlation" %in% igraph::edge_attr_names(g)) {
        edge_strength <- abs(as.numeric(E(g)$correlation))
    } else if ("weight" %in% igraph::edge_attr_names(g)) {
        edge_strength <- abs(as.numeric(E(g)$weight))
    } else {
        edge_strength <- rep(1, igraph::ecount(g))
    }

    if (any(!is.finite(edge_strength))) {
        stop("绘图稀疏化时发现非有限边强度。")
    }

    selected_edge_ids <- integer()

    for (vertex_id in seq_len(igraph::vcount(g))) {
        incident_ids <- as.integer(
            igraph::incident(g, vertex_id, mode = "all")
        )

        if (length(incident_ids) == 0L) {
            next
        }

        # 先按 |r| 降序；强度相同则按原始 edge ID 升序，保证可重复。
        edge_order <- order(
            -edge_strength[incident_ids],
            incident_ids
        )

        keep_number <- min(as.integer(top_k), length(incident_ids))
        selected_edge_ids <- c(
            selected_edge_ids,
            incident_ids[edge_order[seq_len(keep_number)]]
        )
    }

    selected_edge_ids <- sort(unique(selected_edge_ids))

    if (length(selected_edge_ids) == 0L) {
        stop("网页绘图稀疏化后没有保留任何边。")
    }

    plot_g <- igraph::subgraph_from_edges(
        g,
        eids = E(g)[selected_edge_ids],
        delete.vertices = FALSE
    )

    plot_average_degree <- mean(igraph::degree(plot_g))

    list(
        graph = plot_g,
        applied = TRUE,
        original_nodes = igraph::vcount(g),
        original_edges = igraph::ecount(g),
        original_average_degree = original_average_degree,
        plot_nodes = igraph::vcount(plot_g),
        plot_edges = igraph::ecount(plot_g),
        plot_average_degree = plot_average_degree,
        avg_degree_threshold = avg_degree_threshold,
        top_k = as.integer(top_k)
    )
}


# ------------------------------------------------------------------
# 构建可区分正负相关边的 networkD3 交互网络
#
# 说明：
# - 节点颜色仍按指定分类列（默认 Phylum）区分；
# - 正相关边与负相关边分别使用不同颜色；
# - 边宽按 |correlation| 映射；
# - 输入 g 可以是仅供 HTML 展示的稀疏网络副本。
# ------------------------------------------------------------------
build_signed_force_network_widget <- function(
    g,
    node_table,
    node_group_col = "Phylum",
    positive_color = "#D55E00",
    negative_color = "#0072B2"
) {
    if (igraph::vcount(g) == 0L || igraph::ecount(g) == 0L) {
        stop("用于 networkD3 绘图的网络为空。")
    }
    if (!"name" %in% colnames(node_table)) {
        stop("node_table 缺少 name 列。")
    }
    if (!"correlation" %in% igraph::edge_attr_names(g)) {
        stop("绘图网络缺少 correlation 边属性。")
    }

    node_names <- V(g)$name
    node_match <- match(node_names, as.character(node_table$name))
    if (anyNA(node_match)) {
        missing_names <- node_names[is.na(node_match)]
        stop(
            "以下绘图节点未在 node_table 中找到：",
            paste(utils::head(missing_names, 10), collapse = ", ")
        )
    }

    nodes <- node_table[node_match, , drop = FALSE]
    nodes$node_id <- as.character(nodes$name)

    if (node_group_col %in% colnames(nodes)) {
        nodes$node_group <- as.character(nodes[[node_group_col]])
        nodes$node_group[
            is.na(nodes$node_group) | trimws(nodes$node_group) == ""
        ] <- "Unknown"
    } else {
        warning(
            "node_table 中没有 ", node_group_col,
            " 列，networkD3 节点统一归为 Unknown。"
        )
        nodes$node_group <- "Unknown"
    }

    # 保持节点大小基本一致，避免节点大小掩盖正负边颜色。
    nodes$nodesize <- 1

    edge_table <- igraph::as_data_frame(g, what = "edges")
    edge_table$correlation <- as.numeric(edge_table$correlation)
    if (any(!is.finite(edge_table$correlation))) {
        stop("绘图边中存在非有限 correlation。")
    }

    node_index <- setNames(seq_along(node_names) - 1L, node_names)
    links <- data.frame(
        source = unname(node_index[edge_table$from]),
        target = unname(node_index[edge_table$to]),
        value = pmax(abs(edge_table$correlation), 1e-8),
        correlation = edge_table$correlation,
        sign = ifelse(edge_table$correlation > 0, "positive", "negative"),
        link_color = ifelse(
            edge_table$correlation > 0,
            positive_color,
            negative_color
        ),
        stringsAsFactors = FALSE
    )

    widget <- networkD3::forceNetwork(
        Links = links,
        Nodes = nodes,
        Source = "source",
        Target = "target",
        Value = "value",
        NodeID = "node_id",
        Nodesize = "nodesize",
        Group = "node_group",
        zoom = TRUE,
        legend = TRUE,
        linkColour = links$link_color,
        linkWidth = htmlwidgets::JS(
            "function(d) { return 0.7 + 2.3 * Math.sqrt(d.value); }"
        ),
        linkDistance = 50,
        charge = -35,
        fontSize = 7,
        opacity = 0.72,
        opacityNoHover = 0.08
    )

    legend_js <- sprintf(
        paste0(
            "function(el, x) {",
            "var svg = d3.select(el).select('svg');",
            "var boxWidth = 178, boxHeight = 64;",
            "var containerWidth = el.getBoundingClientRect().width;",
            "var xPos = Math.max(containerWidth - boxWidth - 12, 10);",
            "var legend = svg.append('g')",
            ".attr('class', 'signed-edge-legend')",
            ".attr('transform', 'translate(' + xPos + ',12)');",
            "legend.append('rect')",
            ".attr('width', boxWidth).attr('height', boxHeight)",
            ".attr('rx', 4).attr('ry', 4)",
            ".style('fill', 'white').style('fill-opacity', 0.88)",
            ".style('stroke', '#999').style('stroke-width', 0.6);",
            "legend.append('text').attr('x', 10).attr('y', 17)",
            ".style('font-size', '12px').style('font-weight', 'bold')",
            ".text('Edge sign');",
            "legend.append('line').attr('x1', 10).attr('x2', 42)",
            ".attr('y1', 34).attr('y2', 34)",
            ".style('stroke', '%s').style('stroke-width', 3);",
            "legend.append('text').attr('x', 50).attr('y', 38)",
            ".style('font-size', '11px').text('Positive correlation');",
            "legend.append('line').attr('x1', 10).attr('x2', 42)",
            ".attr('y1', 53).attr('y2', 53)",
            ".style('stroke', '%s').style('stroke-width', 3);",
            "legend.append('text').attr('x', 50).attr('y', 57)",
            ".style('font-size', '11px').text('Negative correlation');",
            "}"
        ),
        positive_color,
        negative_color
    )

    htmlwidgets::onRender(widget, legend_js)
}

# ------------------------------------------------------------------
# 将正相关和负相关的 coabundance 图放入同一个 HTML 页面。
# 两个面板仍分别保留各自的矩阵，避免正负连接数相互抵消。
# ------------------------------------------------------------------
write_combined_coabund_html <- function(
    file,
    project,
    positive_html,
    negative_html,
    positive_color = "#D55E00",
    negative_color = "#0072B2"
) {
    positive_name <- basename(positive_html)
    negative_name <- basename(negative_html)

    html <- c(
        "<!doctype html>",
        "<html>",
        "<head>",
        "<meta charset=\"utf-8\">",
        "<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">",
        paste0("<title>", project, " signed coabundance summary</title>"),
        "<style>",
        "*{box-sizing:border-box;}",
        "body{margin:0;padding:14px;background:#f4f5f7;font-family:Arial,Helvetica,sans-serif;color:#111;}",
        ".page{max-width:1500px;margin:0 auto;}",
        "h1{font-size:23px;margin:0 0 6px 0;}",
        ".note{font-size:14px;margin-bottom:14px;color:#444;line-height:1.5;}",
        ".grid{display:grid;grid-template-columns:1fr;gap:18px;}",
        ".panel{background:#fff;border:1px solid #ddd;border-radius:8px;box-shadow:0 1px 3px rgba(0,0,0,.12);overflow:hidden;}",
        paste0(".positive{border-top:5px solid ", positive_color, ";}"),
        paste0(".negative{border-top:5px solid ", negative_color, ";}"),
        ".panel h2{font-size:18px;margin:10px 14px 0 14px;}",
        "iframe{display:block;width:100%;height:1030px;border:0;background:#fff;}",
        "@media(max-width:900px){iframe{height:900px;}body{padding:8px;}}",
        "</style>",
        "</head>",
        "<body>",
        "<main class=\"page\">",
        paste0("<h1>", project, " Phylum-level signed coabundance</h1>"),
        "<div class=\"note\">Positive and negative association counts are shown in separate full-width panels. The two panels use the same Phylum order and colors, while the ribbon outlines distinguish correlation signs.</div>",
        "<div class=\"grid\">",
        "<section class=\"panel positive\">",
        "<h2>Positive correlations</h2>",
        paste0("<iframe src=\"", positive_name, "\" title=\"Positive coabundance\"></iframe>"),
        "</section>",
        "<section class=\"panel negative\">",
        "<h2>Negative correlations</h2>",
        paste0("<iframe src=\"", negative_name, "\" title=\"Negative coabundance\"></iframe>"),
        "</section>",
        "</div>",
        "</main>",
        "</body>",
        "</html>"
    )

    writeLines(html, con = file)
}


# ------------------------------------------------------------------
# 保存 coabundance 方阵。第一列为原始 Phylum 名，便于直接读取和核查。
# ------------------------------------------------------------------
write_coabund_matrix <- function(x, file) {
    if (is.null(x)) {
        write.table(
            data.frame(taxon = character(), stringsAsFactors = FALSE),
            file = file,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE
        )
        return(invisible(NULL))
    }

    out <- data.frame(
        taxon = rownames(x),
        as.data.frame(x, check.names = FALSE),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    write.table(
        out,
        file = file,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE
    )
    invisible(NULL)
}


# ------------------------------------------------------------------
# 对正、负 Phylum 方阵进行定量比较。
# 只取上三角（含对角线），避免对称矩阵中的跨 Phylum 边被重复计数。
# ------------------------------------------------------------------
build_coabund_comparison <- function(positive_matrix, negative_matrix) {
    if (is.null(positive_matrix) || is.null(negative_matrix)) {
        return(list(
            summary = data.frame(
                positive_total_links = NA_real_,
                negative_total_links = NA_real_,
                positive_nonzero_pairs = NA_integer_,
                negative_nonzero_pairs = NA_integer_,
                identical_matrices = NA,
                max_absolute_pair_difference = NA_real_,
                mean_absolute_pair_difference = NA_real_,
                pair_count_correlation = NA_real_,
                nonzero_pair_jaccard = NA_real_,
                stringsAsFactors = FALSE
            ),
            pairs = data.frame()
        ))
    }

    common_taxa <- intersect(rownames(positive_matrix), rownames(negative_matrix))
    positive_matrix <- positive_matrix[common_taxa, common_taxa, drop = FALSE]
    negative_matrix <- negative_matrix[common_taxa, common_taxa, drop = FALSE]

    use_idx <- which(upper.tri(positive_matrix, diag = TRUE), arr.ind = TRUE)
    pos_vals <- as.numeric(positive_matrix[use_idx])
    neg_vals <- as.numeric(negative_matrix[use_idx])

    pos_nonzero <- pos_vals > 0
    neg_nonzero <- neg_vals > 0
    union_nonzero <- pos_nonzero | neg_nonzero
    intersection_nonzero <- pos_nonzero & neg_nonzero

    pair_cor <- NA_real_
    if (length(pos_vals) > 1L && stats::sd(pos_vals) > 0 && stats::sd(neg_vals) > 0) {
        pair_cor <- suppressWarnings(stats::cor(pos_vals, neg_vals, method = "pearson"))
    }

    jaccard <- if (sum(union_nonzero) > 0) {
        sum(intersection_nonzero) / sum(union_nonzero)
    } else {
        NA_real_
    }

    identical_matrices <- isTRUE(all.equal(
        unname(positive_matrix),
        unname(negative_matrix),
        tolerance = 0,
        check.attributes = FALSE
    ))

    pairs <- data.frame(
        taxon_1 = rownames(positive_matrix)[use_idx[, 1]],
        taxon_2 = colnames(positive_matrix)[use_idx[, 2]],
        positive_links = pos_vals,
        negative_links = neg_vals,
        total_links = pos_vals + neg_vals,
        positive_minus_negative = pos_vals - neg_vals,
        absolute_difference = abs(pos_vals - neg_vals),
        stringsAsFactors = FALSE
    )
    pairs <- pairs[order(-pairs$total_links, -pairs$absolute_difference), , drop = FALSE]

    summary <- data.frame(
        positive_total_links = sum(pos_vals, na.rm = TRUE),
        negative_total_links = sum(neg_vals, na.rm = TRUE),
        positive_nonzero_pairs = sum(pos_nonzero, na.rm = TRUE),
        negative_nonzero_pairs = sum(neg_nonzero, na.rm = TRUE),
        identical_matrices = identical_matrices,
        max_absolute_pair_difference = max(abs(pos_vals - neg_vals), na.rm = TRUE),
        mean_absolute_pair_difference = mean(abs(pos_vals - neg_vals), na.rm = TRUE),
        pair_count_correlation = pair_cor,
        nonzero_pair_jaccard = jaccard,
        stringsAsFactors = FALSE
    )

    list(summary = summary, pairs = pairs)
}


# ------------------------------------------------------------------
# 使用同一组 Phylum 和同一套颜色生成正、负相关 chord 图。
# 图中显示名称去除 p__ 前缀，tooltip 与矩阵文件仍保留原始名称。
# 增大 margin、固定画布大小并改为纵向全宽面板，避免长标签被裁切。
# ------------------------------------------------------------------
build_aligned_coabund_widgets <- function(
    positive_matrix,
    negative_matrix,
    plot_num = 10L,
    color_values = RColorBrewer::brewer.pal(10, "Paired"),
    positive_chord_color = "#D55E00",
    negative_chord_color = "#0072B2"
) {
    matrix_taxa <- function(x) {
        if (is.null(x)) {
            return(character())
        }
        unique(c(rownames(x), colnames(x)))
    }

    all_taxa <- sort(unique(c(
        matrix_taxa(positive_matrix),
        matrix_taxa(negative_matrix)
    )))

    if (length(all_taxa) == 0L) {
        return(list(
            positive = NULL,
            negative = NULL,
            taxa = character(),
            display_taxa = character(),
            positive_matrix = NULL,
            negative_matrix = NULL
        ))
    }

    expand_matrix <- function(x, taxa) {
        out <- matrix(
            0,
            nrow = length(taxa),
            ncol = length(taxa),
            dimnames = list(taxa, taxa)
        )
        if (!is.null(x) && nrow(x) > 0L && ncol(x) > 0L) {
            use_rows <- intersect(rownames(x), taxa)
            use_cols <- intersect(colnames(x), taxa)
            out[use_rows, use_cols] <- as.matrix(x[use_rows, use_cols, drop = FALSE])
        }
        out
    }

    positive_full <- expand_matrix(positive_matrix, all_taxa)
    negative_full <- expand_matrix(negative_matrix, all_taxa)

    combined_score <- rowSums(positive_full, na.rm = TRUE) +
        colSums(positive_full, na.rm = TRUE) +
        rowSums(negative_full, na.rm = TRUE) +
        colSums(negative_full, na.rm = TRUE)
    taxa_order <- names(sort(combined_score, decreasing = TRUE))
    taxa_order <- taxa_order[combined_score[taxa_order] > 0]

    if (length(taxa_order) == 0L) {
        return(list(
            positive = NULL,
            negative = NULL,
            taxa = character(),
            display_taxa = character(),
            positive_matrix = NULL,
            negative_matrix = NULL
        ))
    }

    use_number <- min(as.integer(plot_num), length(taxa_order))
    selected_taxa <- taxa_order[seq_len(use_number)]
    display_taxa <- sub("^[A-Za-z]__", "", selected_taxa)

    positive_use <- positive_full[selected_taxa, selected_taxa, drop = FALSE]
    negative_use <- negative_full[selected_taxa, selected_taxa, drop = FALSE]

    if (length(color_values) < length(selected_taxa)) {
        stop(
            "coabundance 颜色数量不足：需要 ", length(selected_taxa),
            " 个，实际只有 ", length(color_values), " 个。"
        )
    }
    selected_colors <- color_values[seq_along(selected_taxa)]

    common_chord_args <- list(
        width = 1200,
        height = 930,
        margin = 210,
        groupNames = display_taxa,
        tooltipNames = selected_taxa,
        groupColors = selected_colors,
        groupPadding = 3,
        groupnamePadding = 18,
        groupnameFontsize = 15,
        ticklabelFontsize = 9,
        showTicks = TRUE,
        showZeroTooltips = FALSE,
        fadeLevel = 0.12
    )

    positive_widget <- NULL
    negative_widget <- NULL

    if (sum(positive_use, na.rm = TRUE) > 0) {
        positive_widget <- do.call(
            chorddiag::chorddiag,
            c(
                list(data = positive_use, chordedgeColor = positive_chord_color),
                common_chord_args
            )
        )
    }

    if (sum(negative_use, na.rm = TRUE) > 0) {
        negative_widget <- do.call(
            chorddiag::chorddiag,
            c(
                list(data = negative_use, chordedgeColor = negative_chord_color),
                common_chord_args
            )
        )
    }

    list(
        positive = positive_widget,
        negative = negative_widget,
        taxa = selected_taxa,
        display_taxa = display_taxa,
        positive_matrix = positive_use,
        negative_matrix = negative_use
    )
}

# ==================== 参数解析 ====================
parser <- OptionParser(description = "Signed Microbial Network Analysis via FastSpar")
parser <- add_option(parser, c("-p", "--project"), help = "Project Name")
parser <- add_option(parser, c("-l", "--level"), help = "Taxonomy Level: species or genus")
parser <- add_option(parser, c("-f", "--filiter"), default = "0.3", help = "Absolute correlation threshold, e.g. 0.3")
parser <- add_option(parser, c("-a", "--output_a"), default = "./", help = "Output Dir A")
parser <- add_option(parser, c("-b", "--output_b"), default = "./", help = "Output Dir B")
parser <- add_option(parser, c("-c", "--output_c"), default = "./", help = "Output Dir C")
parser <- add_option(parser, c("-g", "--output_gexf"), default = "./", help = "Output GEXF Dir")
parser <- add_option(parser, c("-d", "--output_node"), default = "./", help = "Output Node/Edge Table Dir")
parser <- add_option(parser, c("-C", "--cor_matrix"), help = "Path to FastSpar correlation matrix")
parser <- add_option(parser, c("-P", "--p_matrix"), help = "Path to FastSpar p-value matrix")
parser <- add_option(parser, c("--spins"), type = "integer", default = 25L, help = "Maximum number of signed spin-glass communities per connected component [default %default]")
parser <- add_option(parser, c("--gamma"), type = "double", default = 1, help = "Positive-edge resolution parameter [default %default]")
parser <- add_option(parser, c("--gamma-minus"), dest = "gamma_minus", type = "double", default = 1, help = "Negative-edge resolution parameter [default %default]")
parser <- add_option(parser, c("--seed"), type = "integer", default = 12345L, help = "Random seed for signed spin-glass [default %default]")
parser <- add_option(parser, c("--plot-degree-threshold"), dest = "plot_degree_threshold", type = "double", default = 12, help = "If full-network average degree exceeds this value, sparsify only network_plot.html [default %default]")
parser <- add_option(parser, c("--plot-top-k"), dest = "plot_top_k", type = "integer", default = 3L, help = "For HTML visualization, retain each node's top-k strongest |r| edges and take their union [default %default]")
parser <- add_option(parser, c("--positive-edge-color"), dest = "positive_edge_color", type = "character", default = "#D55E00", help = "Positive correlation edge color for network_plot.html [default %default]")
parser <- add_option(parser, c("--negative-edge-color"), dest = "negative_edge_color", type = "character", default = "#0072B2", help = "Negative correlation edge color for network_plot.html [default %default]")
parser <- add_option(parser, c("--species-base"), dest="species_base", default=Sys.getenv("MACKDB_SPECIES_INPUT", "data/species"))
parser <- add_option(parser, c("--genus-base"), dest="genus_base", default=Sys.getenv("MACKDB_GENUS_INPUT", "data/genus"))
parser <- add_option(parser, c("--taxonomy-file"), dest="taxonomy_file", default="data/reference/bac120_taxonomy_r226_species_uniq.tsv")
parse <- parse_args(parser)

required_options <- c("project", "level", "cor_matrix", "p_matrix")
missing_options <- required_options[
    vapply(required_options, function(x) {
        is.null(parse[[x]]) || is.na(parse[[x]]) || parse[[x]] == ""
    }, logical(1))
]
if (length(missing_options) > 0L) {
    stop("缺少必需参数：", paste(missing_options, collapse = ", "))
}
if (parse$spins < 2L) {
    stop("--spins 必须不小于 2。")
}
if (!is.finite(parse$gamma) || parse$gamma < 0) {
    stop("--gamma 必须是非负有限数值。")
}
if (!is.finite(parse$gamma_minus) || parse$gamma_minus < 0) {
    stop("--gamma-minus 必须是非负有限数值。")
}
if (!is.finite(parse$plot_degree_threshold) || parse$plot_degree_threshold < 0) {
    stop("--plot-degree-threshold 必须是非负有限数值。")
}
if (parse$plot_top_k < 1L) {
    stop("--plot-top-k 必须不小于 1。")
}

output_a_dir <- safe_normalize_dir(parse$output_a)
output_b_dir <- safe_normalize_dir(parse$output_b)
output_c_dir <- safe_normalize_dir(parse$output_c)
output_gexf_dir <- safe_normalize_dir(parse$output_gexf)
output_node_dir <- safe_normalize_dir(parse$output_node)

# ==================== 1. 数据读取与字典构建 ====================
if (parse$level == "species") {
    base_dir <- parse$species_base
    otu_file <- file.path(base_dir, "abundance", paste0(parse$project, "_abundance.matrix"))
    proj_tax_file <- file.path(base_dir, "species", paste0(parse$project, "_species.matrix"))
} else if (parse$level == "genus") {
    base_dir <- parse$genus_base
    otu_file <- file.path(base_dir, "abundance", paste0(parse$project, "_abundance.matrix"))
    proj_tax_file <- file.path(base_dir, "genus", paste0(parse$project, "_genus.matrix"))
} else {
    stop("分类层级 -l 必须是 'species' 或 'genus'")
}

otu <- read.csv(otu_file, sep = "\t", check.names = FALSE)
names(otu)[1] <- "key"

proj_tax <- read.csv(
    proj_tax_file,
    sep = "\t",
    col.names = c("key", "replacement"),
    header = FALSE
)

all_tax_species <- read.csv(
    parse$taxonomy_file,
    row.names = 1,
    sep = ";"
)

proj_tax <- proj_tax %>%
    mutate(clean_name = sub("^\\d+_([^_]+.*)$", "\\1", replacement))
id_to_name <- setNames(proj_tax$clean_name, as.character(proj_tax$key))

merged_data <- merge(
    otu,
    proj_tax[, c("key", "clean_name")],
    by = "key",
    all.x = TRUE
)
merged_data <- merged_data[!is.na(merged_data$clean_name), ]
merged_data <- merged_data[!duplicated(merged_data$clean_name), ]
rownames(merged_data) <- merged_data$clean_name
merged_data$clean_name <- NULL
merged_data$key <- NULL

# ==================== 2. 构建 microtable ====================
if (parse$level == "species") {
    tax_sp <- all_tax_species
    clean_sp <- gsub("s__", "", tax_sp$Species)
    clean_sp <- gsub(" ", "_", clean_sp)

    tax_sp$match_name <- clean_sp
    tax_sp <- tax_sp[!duplicated(tax_sp$match_name), ]
    rownames(tax_sp) <- tax_sp$match_name
    tax_sp$match_name <- NULL

    dataset <- microtable$new(
        otu_table = merged_data,
        tax_table = tax_sp
    )
} else {
    all_tax_genus <- all_tax_species[, !(names(all_tax_species) %in% "Species")]
    all_tax_genus <- all_tax_genus[!duplicated(all_tax_genus), ]
    rownames(all_tax_genus) <- gsub("g__", "", all_tax_genus$Genus)

    dataset <- microtable$new(
        otu_table = merged_data,
        tax_table = all_tax_genus
    )
}
dataset$tidy_dataset()

# ==================== 3. 导入 FastSpar 矩阵 ====================
message("正在导入 FastSpar 相关性与 P 值矩阵...")

ext_cor <- read.delim(
    parse$cor_matrix,
    row.names = 1,
    check.names = FALSE
)
ext_p <- read.delim(
    parse$p_matrix,
    row.names = 1,
    check.names = FALSE
)

rownames(ext_cor) <- id_to_name[rownames(ext_cor)]
colnames(ext_cor) <- id_to_name[colnames(ext_cor)]
rownames(ext_p) <- id_to_name[rownames(ext_p)]
colnames(ext_p) <- id_to_name[colnames(ext_p)]

# ==================== 4. 初始化网络并筛选物种 ====================
t1 <- trans_network$new(
    dataset = dataset,
    cor_method = NULL,
    filter_thres = 0.001
)

valid_taxa <- Reduce(
    intersect,
    list(
        colnames(t1$data_abund),
        rownames(ext_cor),
        colnames(ext_cor),
        rownames(ext_p),
        colnames(ext_p)
    )
)

if (length(valid_taxa) < 2L) {
    stop("完成丰度过滤和矩阵对齐后，可用分类单元少于 2 个。")
}

final_cor <- as.matrix(ext_cor[valid_taxa, valid_taxa, drop = FALSE])
final_p <- as.matrix(ext_p[valid_taxa, valid_taxa, drop = FALSE])
storage.mode(final_cor) <- "double"
storage.mode(final_p) <- "double"

# FastSpar 理论上输出对称矩阵。这里在筛边前做数值对称化：
# 相关系数取两侧均值；P 值取两侧较大值（更保守）。
cor_asymmetry <- suppressWarnings(max(abs(final_cor - t(final_cor)), na.rm = TRUE))
p_asymmetry <- suppressWarnings(max(abs(final_p - t(final_p)), na.rm = TRUE))
if (is.finite(cor_asymmetry) && cor_asymmetry > 1e-10) {
    message("相关性矩阵存在轻微不对称，最大差值=", signif(cor_asymmetry, 6), "；将取两侧均值。")
}
if (is.finite(p_asymmetry) && p_asymmetry > 1e-10) {
    message("P值矩阵存在轻微不对称，最大差值=", signif(p_asymmetry, 6), "；将取两侧较大值。")
}
final_cor <- (final_cor + t(final_cor)) / 2
p_left <- final_p
p_right <- t(final_p)
final_p <- ifelse(
    is.na(p_left),
    p_right,
    ifelse(is.na(p_right), p_left, pmax(p_left, p_right))
)
dim(final_p) <- dim(p_left)
dimnames(final_p) <- dimnames(p_left)

# ==================== 5. 根据阈值生成有符号网络 ====================
cor_cut_val <- as.numeric(parse$filiter)
p_thres <- 0.0001

if (!is.finite(cor_cut_val) || cor_cut_val < 0 || cor_cut_val > 1) {
    stop("相关性绝对值阈值 -f 必须位于 [0, 1]。")
}

message(
    "构建有符号网络 -> |相关性|阈值: ", cor_cut_val,
    "；P值阈值: ", p_thres
)

adj_matrix <- final_cor
adj_matrix[abs(adj_matrix) < cor_cut_val] <- 0
adj_matrix[final_p > p_thres] <- 0
diag(adj_matrix) <- 0
adj_matrix[is.na(adj_matrix)] <- 0

# weighted=TRUE 时，初始 E(g)$weight 是带符号相关系数。
g <- igraph::graph_from_adjacency_matrix(
    adj_matrix,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
)

# 永久保存原始相关系数及方向。
E(g)$correlation <- as.numeric(E(g)$weight)
E(g)$sign <- ifelse(E(g)$correlation > 0, "positive", "negative")
E(g)$label <- ifelse(E(g)$correlation > 0, "+", "-")

# weight 仅用于绘图、常规中心性及 microeco 中要求非负权重的函数。
# signed spin-glass 会显式使用 E(g)$correlation，而不是 E(g)$weight。
E(g)$weight <- abs(E(g)$correlation)

# 删除没有任何显著边的节点。
node_num_before <- igraph::vcount(g)
isolated_nodes <- igraph::V(g)[igraph::degree(g) == 0]
if (length(isolated_nodes) > 0L) {
    g <- igraph::delete_vertices(g, isolated_nodes)
}
message("删除孤立节点: ", length(isolated_nodes), " / ", node_num_before)

if (igraph::vcount(g) == 0L || igraph::ecount(g) == 0L) {
    message("阈值筛选并删除孤立节点后无可分析网络，写出空表。")

    empty_node_table <- data.frame(
        name = character(),
        degree = integer(),
        betweenness_centrality = numeric(),
        closeness_centrality = numeric(),
        eigenvector_centrality = numeric(),
        Abundance = numeric(),
        module = character(),
        z = numeric(),
        p = numeric(),
        taxa_roles = character(),
        positive_degree = integer(),
        negative_degree = integer(),
        within_positive_degree = integer(),
        within_negative_degree = integer(),
        between_positive_degree = integer(),
        between_negative_degree = integer(),
        positive_strength = numeric(),
        negative_strength_abs = numeric(),
        stringsAsFactors = FALSE
    )

    empty_edge_table <- data.frame(
        from = character(),
        to = character(),
        weight = numeric(),
        correlation = numeric(),
        sign = character(),
        label = character(),
        from_module = character(),
        to_module = character(),
        same_module = logical(),
        stringsAsFactors = FALSE
    )

    write.table(
        empty_node_table,
        file = file.path(output_node_dir, paste0(parse$project, "_node.txt")),
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE
    )
    write.table(
        empty_edge_table,
        file = file.path(output_node_dir, paste0(parse$project, "_edge_signed.txt")),
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE
    )

    message("网络构建与保存完成！")
    quit(save = "no", status = 0)
}

# 手动注入丰度信息供后续绘图使用。
node_abund_mat <- t1$data_abund[, V(g)$name, drop = FALSE]
mean_abund <- apply(node_abund_mat, 2, mean)
V(g)$Abundance <- mean_abund[V(g)$name]
V(g)$RelativeAbundance <- V(g)$Abundance / sum(V(g)$Abundance)

# ==================== 6. 有符号模块划分 ====================
message(
    "开始有符号模块划分：implementation=neg；spins=", parse$spins,
    "；gamma=", parse$gamma,
    "；gamma.minus=", parse$gamma_minus,
    "；seed=", parse$seed
)

signed_module_result <- assign_signed_spinglass_modules(
    g = g,
    spins = parse$spins,
    gamma = parse$gamma,
    gamma_minus = parse$gamma_minus,
    seed = parse$seed
)
g <- signed_module_result$graph
component_summary <- signed_module_result$component_summary

message("有符号模块划分完成，共识别模块数: ", length(unique(V(g)$module)))

# 将已带有 signed module 的网络注入 microeco。
t1$res_network <- g

# 不调用 t1$cal_module()，防止用普通社区算法覆盖 signed module。
# 不调用 t1$cal_network_attr()，因为其内部固定使用 cluster_fast_greedy，
# 返回的是基于 abs(weight) 的普通 modularity，并非本次 signed module 的质量指标。

# ==================== 7. 标准 Zi/Pi 与节点属性 ====================
message("根据 signed module 计算标准 Zi/Pi 和节点属性...")
node_roles_ok <- TRUE

tryCatch(
    {
        t1$get_node_table(node_roles = TRUE)
        cs <- t1$res_node_table
    },
    error = function(e) {
        node_roles_ok <<- FALSE
        message("节点角色计算失败，改为写出基础节点表：", conditionMessage(e))
        cs <<- build_basic_node_table(g, dataset)
        t1$res_node_table <<- cs
    }
)

signed_node_metrics <- calculate_signed_node_metrics(g)
cs <- dplyr::left_join(cs, signed_node_metrics, by = "name")
rownames(cs) <- cs$name
t1$res_node_table <- cs

# 输出边与模块诊断表。
signed_edge_table <- build_signed_edge_table(g)
module_summary <- build_module_summary(g, signed_edge_table)

message(
    "最终网络：节点=", igraph::vcount(g),
    "；边=", igraph::ecount(g),
    "；正边=", sum(E(g)$correlation > 0),
    "；负边=", sum(E(g)$correlation < 0),
    "；连通分量=", igraph::components(g)$no,
    "；模块=", length(unique(V(g)$module))
)

# ==================== 8. 优先保存核心文本结果 ====================
# 核心表格先于绘图和 GEXF 保存。这样即使后续 HTML/GEXF 导出失败，
# node、edge、module 和 component 结果也不会丢失。
write.table(
    cs,
    file = file.path(output_node_dir, paste0(parse$project, "_node.txt")),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
)

write.table(
    signed_edge_table,
    file = file.path(output_node_dir, paste0(parse$project, "_edge_signed.txt")),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
)

write.table(
    module_summary,
    file = file.path(output_node_dir, paste0(parse$project, "_module_signed_summary.txt")),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
)

write.table(
    component_summary,
    file = file.path(output_node_dir, paste0(parse$project, "_component_signed_summary.txt")),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
)

message("核心节点、边和模块表已保存。")

# ==================== 9. 构建仅供 HTML 展示的稀疏网络 ====================
plot_network_info <- build_plot_network(
    g = g,
    avg_degree_threshold = parse$plot_degree_threshold,
    top_k = parse$plot_top_k
)
g_plot <- plot_network_info$graph

if (isTRUE(plot_network_info$applied)) {
    message(
        "network_plot.html 展示稀疏化已启用：完整网络平均度=",
        round(plot_network_info$original_average_degree, 4),
        " > ", plot_network_info$avg_degree_threshold,
        "；每个节点保留 |r| 最强的 ", plot_network_info$top_k,
        " 条边并取并集；边数 ", plot_network_info$original_edges,
        " -> ", plot_network_info$plot_edges,
        "；展示网络平均度=", round(plot_network_info$plot_average_degree, 4),
        "。完整网络仍用于 signed module、Zi/Pi、文本结果和 GEXF。"
    )
} else {
    message(
        "network_plot.html 使用完整网络：平均度=",
        round(plot_network_info$original_average_degree, 4),
        " <= 阈值 ", plot_network_info$avg_degree_threshold,
        "。"
    )
}

plot_network_summary <- data.frame(
    project = parse$project,
    strategy = ifelse(
        plot_network_info$applied,
        "top_k_absolute_correlation_union",
        "full_network"
    ),
    only_for_html_visualization = TRUE,
    average_degree_threshold = plot_network_info$avg_degree_threshold,
    top_k_per_node = plot_network_info$top_k,
    original_nodes = plot_network_info$original_nodes,
    original_edges = plot_network_info$original_edges,
    original_average_degree = plot_network_info$original_average_degree,
    plot_nodes = plot_network_info$plot_nodes,
    plot_edges = plot_network_info$plot_edges,
    plot_positive_edges = sum(E(g_plot)$correlation > 0),
    plot_negative_edges = sum(E(g_plot)$correlation < 0),
    plot_average_degree = plot_network_info$plot_average_degree,
    positive_edge_color = parse$positive_edge_color,
    negative_edge_color = parse$negative_edge_color,
    stringsAsFactors = FALSE
)

write.table(
    plot_network_summary,
    file = file.path(
        output_node_dir,
        paste0(parse$project, "_network_plot_summary.txt")
    ),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
)

# ==================== 10. 图表绘制 ====================
if (node_roles_ok && nrow(cs) > 0L && any(!is.na(cs$z) & !is.na(cs$p))) {
    p_1 <- tryCatch(
        t1$plot_taxa_roles(
            roles_color_background = TRUE,
            add_label = TRUE,
            add_label_group = "Network hubs",
            add_label_text = "name",
            label_text_size = 4,
            label_text_color = "grey50",
            label_text_parse = TRUE,
            plot_module = TRUE,
            use_type = 1
        ),
        error = function(e) {
            message("TopoMap use_type=1 绘图失败：", conditionMessage(e))
            empty_plot("TopoMap unavailable")
        }
    )

    p_2 <- tryCatch(
        t1$plot_taxa_roles(
            use_type = 2,
            use_level = "Phylum",
            show_value = c("z", "p"),
            show_number = 1:10,
            plot_color = "Phylum",
            plot_shape = "taxa_roles",
            plot_size = "Abundance",
            color_values = RColorBrewer::brewer.pal(12, "Paired")
        ),
        error = function(e) {
            message("TopoMap use_type=2 绘图失败：", conditionMessage(e))
            empty_plot("Taxa roles unavailable")
        }
    )
} else {
    p_1 <- empty_plot("No valid Zi/Pi roles")
    p_2 <- empty_plot("No valid Zi/Pi roles")
}

# network_plot.html 使用稀疏网络副本，但正负边符号和原始相关系数均保留。
p1 <- tryCatch(
    build_signed_force_network_widget(
        g = g_plot,
        node_table = cs,
        node_group_col = "Phylum",
        positive_color = parse$positive_edge_color,
        negative_color = parse$negative_edge_color
    ),
    error = function(e) {
        message("signed networkD3 绘图失败：", conditionMessage(e))
        NULL
    }
)

# 共丰度汇总始终使用完整网络，而不是 HTML 稀疏网络。
t1$res_network <- g
t1$res_edge_table <- NULL

coabund_ok <- tryCatch(
    {
        t1$cal_sum_links(taxa_level = "Phylum")
        TRUE
    },
    error = function(e) {
        message("coabundance 汇总计算失败：", conditionMessage(e))
        FALSE
    }
)

p2_positive <- NULL
p2_negative <- NULL
coabund_selected_taxa <- character()

if (isTRUE(coabund_ok)) {
    signed_coabund_widgets <- tryCatch(
        build_aligned_coabund_widgets(
            positive_matrix = t1$res_sum_links_pos,
            negative_matrix = t1$res_sum_links_neg,
            plot_num = 10L,
            color_values = RColorBrewer::brewer.pal(10, "Paired"),
            positive_chord_color = parse$positive_edge_color,
            negative_chord_color = parse$negative_edge_color
        ),
        error = function(e) {
            message("正负 coabund 对齐绘图失败：", conditionMessage(e))
            NULL
        }
    )

    if (!is.null(signed_coabund_widgets)) {
        p2_positive <- signed_coabund_widgets$positive
        p2_negative <- signed_coabund_widgets$negative
        coabund_selected_taxa <- signed_coabund_widgets$taxa

        message(
            "正负 coabund 使用相同的 Phylum 集合：",
            paste(coabund_selected_taxa, collapse = ", ")
        )

        # 保存实际用于绘图的正、负 Phylum 方阵，并定量判断两者是否相同。
        write_coabund_matrix(
            signed_coabund_widgets$positive_matrix,
            file.path(output_node_dir, paste0(parse$project, "_coabund_positive_matrix.txt"))
        )
        write_coabund_matrix(
            signed_coabund_widgets$negative_matrix,
            file.path(output_node_dir, paste0(parse$project, "_coabund_negative_matrix.txt"))
        )

        coabund_comparison <- build_coabund_comparison(
            signed_coabund_widgets$positive_matrix,
            signed_coabund_widgets$negative_matrix
        )

        write.table(
            coabund_comparison$summary,
            file = file.path(output_node_dir, paste0(parse$project, "_coabund_comparison_summary.txt")),
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE
        )
        write.table(
            coabund_comparison$pairs,
            file = file.path(output_node_dir, paste0(parse$project, "_coabund_signed_pair_summary.txt")),
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE
        )

        message(
            "coabund 正负矩阵核查：positive_total=",
            coabund_comparison$summary$positive_total_links,
            "；negative_total=",
            coabund_comparison$summary$negative_total_links,
            "；identical=",
            coabund_comparison$summary$identical_matrices,
            "；pair_count_correlation=",
            round(coabund_comparison$summary$pair_count_correlation, 4)
        )
    }
}

# ==================== 11. 图形与网络文件保存 ====================
plot_combined <- arrangeGrob(p_1, p_2, ncol = 2)
ggsave(
    file = file.path(output_a_dir, paste0(parse$project, "_TopoMap.png")),
    plot = plot_combined,
    width = 20,
    height = 9
)

# R 启动时需要 conda 的旧动态库；pandoc 保存 self-contained html 时不能继承它。
Sys.unsetenv("LD_LIBRARY_PATH")

network_html <- file.path(
    output_b_dir,
    paste0(parse$project, "_network_plot.html")
)
coabund_positive_html <- file.path(
    output_c_dir,
    paste0(parse$project, "_coabund_positive.html")
)
coabund_negative_html <- file.path(
    output_c_dir,
    paste0(parse$project, "_coabund_negative.html")
)
coabund_combined_html <- file.path(
    output_c_dir,
    paste0(parse$project, "_coabund.html")
)

if (is.null(p1)) {
    write_empty_html(network_html, "signed networkD3 plot unavailable")
} else {
    saveWidget(p1, file = network_html, selfcontained = TRUE)
}

if (is.null(p2_positive)) {
    write_empty_html(
        coabund_positive_html,
        "positive coabundance summary plot unavailable"
    )
} else {
    saveWidget(
        p2_positive,
        file = coabund_positive_html,
        selfcontained = TRUE
    )
}

if (is.null(p2_negative)) {
    write_empty_html(
        coabund_negative_html,
        "negative coabundance summary plot unavailable"
    )
} else {
    saveWidget(
        p2_negative,
        file = coabund_negative_html,
        selfcontained = TRUE
    )
}

write_combined_coabund_html(
    file = coabund_combined_html,
    project = parse$project,
    positive_html = coabund_positive_html,
    negative_html = coabund_negative_html,
    positive_color = parse$positive_edge_color,
    negative_color = parse$negative_edge_color
)

message(
    "coabundance HTML 已保存：正相关=", coabund_positive_html,
    "；负相关=", coabund_negative_html,
    "；正负双面板=", coabund_combined_html
)

# GEXF 中：
# - weight 是 abs(correlation)，便于 Gephi 使用正权重布局/边宽；
# - correlation 保留原始正负值；
# - sign 与 label 保留方向；
# - module 是 signed spin-glass 模块。
#
# microeco 2.0.0 的 save_network() 会把 igraph 的所有 graph-level
# attributes 追加到 rgexf 的 meta 列表；而 rgexf 的 meta 只接受
# creator、description、keywords 三个名称。因此导出临时副本时删除
# graph-level attributes，但不删除任何节点或边属性。
gexf_file <- file.path(
    output_gexf_dir,
    paste0(parse$project, ".gexf")
)

g_for_gexf <- g
graph_level_attributes <- igraph::graph_attr_names(g_for_gexf)
if (length(graph_level_attributes) > 0L) {
    message(
        "GEXF 导出前从临时网络副本移除 graph-level attributes：",
        paste(graph_level_attributes, collapse = ", ")
    )
    for (attr_name in graph_level_attributes) {
        g_for_gexf <- igraph::delete_graph_attr(g_for_gexf, attr_name)
    }
}

original_network <- t1$res_network
t1$res_network <- g_for_gexf

gexf_ok <- tryCatch(
    {
        t1$save_network(filepath = gexf_file)
        TRUE
    },
    error = function(e) {
        message("GEXF 导出失败，但核心文本结果已保存：", conditionMessage(e))
        FALSE
    },
    finally = {
        # 恢复包含 signed module 参数的原始网络对象。
        t1$res_network <- original_network
    }
)

if (isTRUE(gexf_ok)) {
    message("GEXF 文件已保存：", gexf_file)
} else {
    # GEXF 失败时额外保存 GraphML；Gephi 和 Cytoscape 均可读取。
    graphml_file <- file.path(
        output_gexf_dir,
        paste0(parse$project, ".graphml")
    )
    tryCatch(
        {
            igraph::write_graph(g, graphml_file, format = "graphml")
            message("已改为保存 GraphML：", graphml_file)
        },
        error = function(e) {
            message("GraphML 备用导出也失败：", conditionMessage(e))
        }
    )
}

message("有符号网络构建、模块划分、Zi/Pi计算与保存完成！")
