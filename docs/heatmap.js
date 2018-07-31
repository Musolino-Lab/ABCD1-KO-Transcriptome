

let range = n => [...Array(n).keys()];
let clamp = (min, max) => ((x) => Math.min(Math.max(x, min), max));
let transpose = (array) => array[0].map((col, i) => array.map(row => row[i]));
let flatten = (array) => [].concat.apply([], array);

function Heatmap(samples_by_genes_matrix, gene_sets, classes) {

    // samples_by_genes_matrix: {'sample_id': {'gene1': float, 'gene2': float}};

    // gene_sets: {'gene_set_name': ['gene1', 'gene2']};

    // classes: {'sample_id': {'category1': 'value', 'category2': value}}

    // categories: ['column_name_of_sample_class_labels_1', 'column_name_of_sample_class_labels_2']
    var categories = _.uniq(_.flatten(Object.values(classes).map(obj => Object.keys(obj))));

    // categories_to_members_to_values: {'category': {'sample_id1': 'value', 'sample_id2': value}}
    var categories_to_members_to_values = _.object(categories.map((category) =>
        [category, _.object(Object.entries(classes).map(([sample_id, categories_to_values]) => [sample_id, categories_to_values[category]]))]
    ));

    // categories_to_values_to_members: {'category': {'value1': ['sample_id1', 'sample_id2']}}  <- I need this one
    var categories_to_values_to_members = _(categories_to_members_to_values).mapObject((samples_to_values) =>
        _(_(Object.entries(samples_to_values)).groupBy(([sample_id, value]) => value)).mapObject((arr_of_arr) => arr_of_arr.map(arr => arr[0]))
    );

    /////////////////////////////////////////////////////////////////////////////
                    ///////    Structure Variables    ///////
    /////////////////////////////////////////////////////////////////////////////

    var margin = {top: 100, right: 100, bottom: 0, left: 100};

    var w = window.innerWidth - (margin.left + margin.right);
    var h = window.innerHeight - (margin.top + margin.bottom);

    var values = 'counts';
    var show_averages = false;
    var scaling = 'mean';
    var sorting = 'complete';
    var reordering = true;
    var minimum_nonzero = 0;

    var scale_types = {
        mean     : () => d3.scaleLinear(),
        absolute : () => d3.scaleLinear(),
        log      : () => d3.scaleLog().clamp(true),  // in the future this will be a magic scale: https://stackoverflow.com/questions/51584888/how-to-write-a-d3-positive-and-negative-logarithmic-scale
    };
    var domains = {};
    var ranges = {};

    value_accessors = {
        counts: (gene) => gene.counts,
        zscore_stddev: (gene) => gene.samples.map((sample) => sample.zscore_stddev),
        zscore_mad: (gene) => gene.samples.map((sample) => sample.zscore_mad),
    }
    value_accessor = value_accessors.counts;


    /////////////////////////////////////////////////////////////////////////////
                    ///////    Styling Variables    ///////
    /////////////////////////////////////////////////////////////////////////////

    var negative_color = '#0000cc';
    var middle_color = '#c0c0c0';
    var positive_color = '#cc0000';
    var show_legends = false;


    var colors20 = ["#3366cc", "#dc3912", "#ff9900", "#109618", "#990099", "#0099c6", "#dd4477",
                    "#66aa00", "#b82e2e", "#316395", "#994499", "#22aa99", "#aaaa11", "#6633cc",
                    "#e67300", "#8b0707", "#651067", "#329262", "#5574a6", "#3b3eac"];

    var id_colors = d3.scaleOrdinal(d3.schemeCategory10);
    var class_colors = _.object(categories, categories.map(c => d3.scaleOrdinal(d3.schemeCategory10)));
    let color_range = (gene, attr) => {
        samples = _.object(_(gene_wise).findWhere({'gene':gene}).samples.map(sample => [sample.sample, sample[attr]]));
        domain = d3.extent(Object.values(samples)); domain.splice(1, 0, 0);
        scale = d3.scaleLinear().domain(domain).range([negative_color, middle_color, positive_color]);
        return [(sample) => scale(samples[sample]), scale]
    }
    var zscore_stddev_colors = null;
    var zscore_stddev_color_scale = null;
    var zscore_mad_colors = null;
    var zscore_mad_color_scale = null;
    var samples_pc1_domain = [];
    var pc1_colors = null;

    var rect_colors = {
        zscore_stddev: (d) => zscore_stddev_colors(d.sample),
        zscore_mad: (d) => zscore_mad_colors(d.sample),
        pc1: (d) => pc1_colors(d.pc1),
    }

    var spacing = 2;
    var rect_side_length = 16;
    var s = () => rect_side_length+spacing;

    /////////////////////////////////////////////////////////////////////////////
                          ///////    Set Up Chart    ///////
    /////////////////////////////////////////////////////////////////////////////

    var svg = d3.select("#graph-container").append("svg").attr("xmlns", "http://www.w3.org/2000/svg").attr("xmlns:xlink", "http://www.w3.org/1999/xlink");
    var g = svg.append("g").attr("transform", "translate("+margin.left+","+margin.top+")");
    svg.style("cursor", "move");

    var gene_set_name = "";
    var genes = [];
    var gene_wise = [];
    var gene_wise_indexer = {};
    var ordered_gene_wise = [];
    var sample_wise = [];
    var x_scales = {};
    var genes_pc1 = [];
    var samples_pc1 = [];

    var x_scales = {};
    var y = (index) => axis_spacing*index;

    function idFunc(d) { return d ? d.id : this.id; }

    var axes = g.selectAll(".axis");
    var dots = g.selectAll(".dots");
    var line = g.selectAll(".line");
    var halo = g.selectAll(".halo");

    var title = g.append("text")
                 .attr("class", "title")
                 .attr("font-family", "sans-serif")
                 .text("")
                 .style("text-anchor", "middle")
                 .attr("x", 0)
                 .attr("y", 0)
                 .attr("dy", "-3em");

    var legend_lines = g.append("g").attr("class", "legend legendLines");
    var legend_points = g.append("g").attr("class", "legend legendPoints");


    /////////////////////////////////////////////////////////////////////////////
                          ///////    Methods    ///////
    /////////////////////////////////////////////////////////////////////////////

    function restart({selected_gene_set_name=gene_set_name,
                      selected_genes=genes}={}) {

        gene_set_name = selected_gene_set_name;
        genes = selected_genes;
        matrix = _(samples_by_genes_matrix).mapObject((sample) => _(sample).pick(genes));

        sample_wise = Object.entries(matrix).map(([sample, genes]) =>
            Object.entries(genes).map(([gene, count]) => { return {
                'id'     : sample+"_"+gene,
                'sample' : sample,
                'gene'   : gene,
                'count'  : count,
            }})
         );

        gene_wise = transpose(sample_wise).map((by_gene) => {
            min    = d3.min(by_gene, gene => gene.count);
            max    = d3.max(by_gene, gene => gene.count);
            mean   = d3.mean(by_gene, gene => gene.count);
            stddev = d3.deviation(by_gene, gene => gene.count);
            mad    = d3.mean(by_gene.map(gene => Math.abs(mean - gene.count)));
            return by_gene.map((gene) => Object.assign(gene, {
                'min'      : min,
                'max'      : max,
                'mean'     : mean,                  // should we only be counting non-zeros?
                'stddev'   : stddev,
                'mad'      : mad,
                'zscore_stddev' : stddev === 0 ? 0 : (gene.count - mean) / stddev,
                'zscore_mad'    : mad === 0 ? 0    : (gene.count - mean) / mad,
            }));
        });

        gene_wise_indexer =  _.object(gene_wise.map((gene, i) => [gene[0].gene, i]));
        sample_wise_indexer = _.object(sample_wise.map((sample, i) => [sample[0].sample, i]));

        ordered_gene_wise = gene_wise;

        // order();
        render();

    }

    function order({values_=values,
                    sorting_=sorting,
                    minimum_nonzero_=minimum_nonzero}={}) {

        values = values_;
        value_accessor = value_accessors[values];
        sorting = sorting_;
        minimum_nonzero = minimum_nonzero_;
        if (minimum_nonzero > 0) { gene_set_name = ""; }

        // THIS FUNCTION NOW ONLY NEEDS TO PRODUCE ARRAYS:
        // ordered_sample_ids
        // ordered_gene_ids

        // Filter by number non-zeros
        ordered_gene_wise = gene_wise.filter((gene) => gene.num_nonzeros >= minimum_nonzero);

        if (ordered_gene_wise.length === 0) { sample_wise = []; return }

        // Set Gene-Wise PC1
        genes_pc1 = reorder.pca1d(reorder.transpose(ordered_gene_wise.map(value_accessor)));
        _(_.zip(ordered_gene_wise, genes_pc1)).each((gene, pc1) => gene.pc1 = pc1);

        // Order Genes
        if (reordering && ordered_gene_wise.length > 1) {

            if (sorting === 'pc1') {
                permutation_order = reorder.sort_order(genes_pc1);
            } else if (sorting === 'complete') {
                permutation_order = reorder.optimal_leaf_order()(ordered_gene_wise.map(value_accessor)).reverse();  // get dendogram out?
            } else { console.log(' this should never happen '); }

            ordered_gene_wise = reorder.stablepermute(ordered_gene_wise, permutation_order);

        } else { permutation_order = range(ordered_gene_wise.length); }

        // Build Sample-Wise
        sample_wise = Object.entries(matrix).map(([sample, genes]) => { return {
                        'id'     : sample,
                        'sample' : sample,
                        'classes': classes[sample],
                        'pc1'    : null,
                        'genes'  : reorder.stablepermute(
                                    Object.entries(genes).map(([gene, count]) => { return {
                                      'id'       : gene+"_"+sample,
                                      'gene'     : gene,
                                      'value'    : values === 'counts' ? count : _(gene_wise[gene_wise_indexer[gene]].samples).findWhere({'sample':sample})[values],
                                      'num_nonzeros': gene_wise[gene_wise_indexer[gene]].num_nonzeros,
                                    }}).filter((gene) => gene.num_nonzeros >= minimum_nonzero),
                                    permutation_order)
                      }});

        // Set Sample-Wise PC1
        samples_pc1 = reorder.pca1d(ordered_gene_wise.map(value_accessor));
        _(_.zip(sample_wise, samples_pc1)).each(([sample, pc1]) => { sample.pc1 = pc1; });
        samples_pc1_domain = d3.extent(samples_pc1);  // hack TODO get rid of this.

        // Add Averages
        if (show_averages) {
            console.log('show averages! append to sample wise! ');
            // sample_wise.append()
        }

        // Domains and Ranges
        recompute_domains_and_ranges();

    }

    function render({spacing_=spacing}={}) {

        spacing = spacing_;

        title.text(gene_set_name);

        x = d3.scalePoint().domain(ordered_sample_ids).range([0,ordered_sample_ids.length*s]);
        y = d3.scalePoint().domain(ordered_gene_ids).range([0,ordered_gene_ids.length*s]);

        rect = g.selectAll(".rect").data(flatten(gene_wise), idFunc);
        gSym = g.selectAll(".gSym").data(ordered_gene_ids, (d) => d);
        sNam = g.selectAll(".sNam").data(ordered_sample_ids, (d) => d);
        // gGru = g.selectAll(".gGru").data(gene_groups, (d) => d[0]);
        // meta = g.selectAll(".meta").data(Object.entries(categories_to_members_to_values), (d) => d[0]);

        // phase 1
            // rectangles which are exiting fade out
            // gene names which are exiting fade out
            // gene groups which are exiting fade out
        t_last = d3.transition().duration(200);
        if (rect.exit().size() > 0) {  // if lines are already drawn,
            rect.exit().transition(t_last).style("opacity", 0).remove();
            gSym.exit().transition(t_last).style("opacity", 0).remove();

            t_last = t_last.transition().duration(500);
        }
        // phase 2
            // rectangles which are staying get re-arranged
            // gene names which are staying get re-arranged
            // gene groups which are staying get re-arranged
        rect

        // phase 3
            // rectangles which are entering get appended
            // gene names which are entering get appended
            // gene groups which are entering get appended
        rect.enter()


    }

    function style({negative_color_=negative_color,
                    middle_color_=middle_color,
                    positive_color_=positive_color,
                    show_legends_=show_legends}={}) {

        negative_color = negative_color_,
        middle_color = middle_color_,
        positive_color = positive_color_,
        show_legends=show_legends_;


        if (line_coloring_system === 'zscore_stddev' || line_coloring_system === 'zscore_mad') {

            g.selectAll(".text").style("font-weight", 300).on("click", (d) => style({lines_color_by_: d.id}));
            if (!_(ordered_gene_wise).findWhere({'gene': lines_color_by})) { lines_color_by = ordered_gene_wise[0].id }
            d3.select('.text#'+lines_color_by).style("font-weight", 700);

            zscore_stddev_colors = color_range(lines_color_by, 'zscore_stddev');
            zscore_stddev_color_scale = zscore_stddev_colors[1]; zscore_stddev_colors = zscore_stddev_colors[0];
            zscore_mad_colors = color_range(lines_color_by, 'zscore_mad');
            zscore_mad_color_scale = zscore_mad_colors[1]; zscore_mad_colors = zscore_mad_colors[0];

        } else if (line_coloring_system === 'class') {

            lines_color_by = categories[0];  // TODO FIX THIS SHIT

        } else if (line_coloring_system === 'pc1') {

            pc1_colors = d3.scaleLinear().domain(samples_pc1_domain).range([negative_color, positive_color]);
        }

        // dots are bound to ordered_gene_wise.samples
        g.selectAll(".rect")
            .attr("visibility", show_rects ? "visible" : "hidden")
            .style("fill", (d) => rect_colors[rect_coloring_system](d))
            .style("opacity", rect_opacity);

        if (show_legends) {


        } else {

        }

    }


    /////////////////////////////////////////////////////////////////////////////
                          ///////    Drag Axes    ///////
    /////////////////////////////////////////////////////////////////////////////

    function drag_axis_start(d) {
        g.selectAll(".line,.halo").transition().duration(100).style('stroke-opacity', 0);
    }

    function drag_axis(d) {
        dragged_index = _(ordered_gene_wise).findIndex((gene) => gene.id === d[0]);

        let expr = (current_index) => {
            if (current_index < dragged_index) {
                if (current_index < ((d3.event.y-axis_spacing/2) / axis_spacing)) {
                        return y(current_index);
                    } else {
                        return y(current_index) + axis_spacing;
                    }
            } else {  // current_index >= dragged_index
                if (current_index < ((d3.event.y+axis_spacing/2) / axis_spacing)) {
                        return y(current_index) - axis_spacing;
                    } else {
                        return y(current_index);
                    }
            }
        }

        g.selectAll(".axis").attr("transform", function(d, i) { return "translate(0," +  expr(i) + ")" });
        g.selectAll(".dots").attr("transform", function(d, i) { return "translate(0," + (expr(i) - max_point_radius) + ")" });

        d3.select(this).attr("transform", "translate(0," + d3.event.y + ")");
        d3.select(".dots#"+d[0]).attr("transform", (d, i) => "translate(0," + (d3.event.y - max_point_radius) + ")");

    }

    function drag_axis_end(d) {

        dragged_index = _(ordered_gene_wise).findIndex((gene) => gene.id === d[0]);
        old_index = dragged_index;
        new_index = clamp(0, ordered_gene_wise.length)(Math.round(d3.event.y / axis_spacing));

        ordered_gene_wise.splice(new_index, 0, ordered_gene_wise.splice(old_index, 1)[0]);
        sample_wise.forEach((sample) => sample.genes.splice(new_index, 0, sample.genes.splice(old_index, 1)[0]));

        render();

    }

    /////////////////////////////////////////////////////////////////////////////
                          ///////    Brush Axes    ///////
    /////////////////////////////////////////////////////////////////////////////


    function setFocus(d) {
    }

    function removeFocus(d) {
    }

    function GeneCards(d) { window.open("http://www.genecards.org/cgi-bin/carddisp.pl?gene="+d.id,'_blank') }


    /////////////////////////////////////////////////////////////////////////////
                          ///////   Zoom & Resize    ///////
    /////////////////////////////////////////////////////////////////////////////

    svg.call(d3.zoom()
               .scaleExtent([1 / 8, 8])
               .on("zoom", zoomed));

    function zoomed() { g.attr("transform", d3.event.transform); }

    function resize() {
        svg.attr("width", window.innerWidth).attr("height", window.innerHeight);
        w = window.innerWidth - (margin.left + margin.right);
        h = window.innerHeight - (margin.top + margin.bottom);
    }

    d3.select(window).on("resize", resize)

    resize();


    /////////////////////////////////////////////////////////////////////////////
                          ///////    Return    ///////
    /////////////////////////////////////////////////////////////////////////////

    return {
        'restart'     : restart,
        'render'      : render,
        'order'       : order,
        'style'       : style,

        get_sorted_gene_list: () => ordered_gene_wise.length ? _(ordered_gene_wise).pluck('gene') : _(gene_wise).pluck('gene'),

        set_reordering: (reordering_) => { reordering = reordering_; if (reordering) { order(); render(); } },
    }

}





