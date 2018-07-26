

function Braid(samples_by_genes_matrix, gene_sets, classes) {

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
                          ///////    Variables    ///////
    /////////////////////////////////////////////////////////////////////////////

    var margin = {top: 100, right: 100, bottom: 0, left: 100};

    var w = window.innerWidth - (margin.left + margin.right);
    var h = window.innerHeight - (margin.top + margin.bottom);

    var values = 'counts';
    var sorting = 'complete';
    var show_all_zero = true;

    var show_points = true;
    var point_coloring_system = 'identity';
    var points_color_by = categories[0];
    var default_point_color = '#333333';
    var show_lines = true;
    var show_halos = true;
    var line_coloring_system = 'identity';
    var lines_color_by = categories[0];
    var default_line_color = '#333333';
    var negative_color = '#0000cc';
    var middle_color = '#c0c0c0';
    var positive_color = '#cc0000';
    var curve = 'linear';
    var halo_opacity = 0.2;


    var absolute_axis = true;
    var relative_axis = false;
    var margin_between_axes = 10;

    var line_bendiness = 0.5;

    var point_radius = 3;
    var gene_spacing = 50;



    var colors20 = ["#3366cc", "#dc3912", "#ff9900", "#109618", "#990099", "#0099c6", "#dd4477",
                    "#66aa00", "#b82e2e", "#316395", "#994499", "#22aa99", "#aaaa11", "#6633cc",
                    "#e67300", "#8b0707", "#651067", "#329262", "#5574a6", "#3b3eac"];

    var scalings = {
        max: (obj) => obj.max * 1.5,
        mean: (obj) => obj.mean * 2,
    };

    var curves = {
        'linear'   : d3.curveLinear,
        'natural'  : d3.curveNatural,
        'catmull'  : d3.curveCatmullRom.alpha(line_bendiness),
        'monotone' : d3.curveMonotoneY,
        'cardinal' : d3.curveCardinal.tension(1-line_bendiness),
    }


    var id_colors = d3.scaleOrdinal(d3.schemeCategory10);
    var class_colors = _.object(categories, categories.map(c => d3.scaleOrdinal(d3.schemeCategory10)));
    let color_range = (gene, attr) => {
        samples = _.object(_(gene_wise).findWhere({'gene':gene}).samples.map(sample => [sample.sample, sample[attr]]));
        domain = d3.extent(Object.values(samples)); domain.splice(1, 0, 0);
        scale = d3.scaleLinear().domain(domain).range([negative_color, middle_color, positive_color]);
        return (sample) => scale(samples[sample])
    }
    var samples_pc1_domain = [];
    var pc1_colors = null;

    var point_colors = {
        identity: (d) => id_colors(d.sample),
        class: (d) => class_colors[points_color_by](d.classes[points_color_by]),
        unform: (d) => default_point_color,
    }

    var line_colors = {
        identity: (d) => id_colors(d.sample),
        class: (d) => class_colors[lines_color_by](d.classes[lines_color_by]),
        unform: (d) => default_line_color,
        z_stddev: (d) => color_range(lines_color_by, 'zscore')(d.sample),
        z_mad: (d) => color_range(lines_color_by, 'zscore_mad')(d.sample),
        pc1: (d) => pc1_colors(d.pc1),
    }

    var line_from_sample = d3.line()
                              .y((d, i) => y(i))
                              .x((d) => x_scales[d.gene](d.count))
                              .curve(curves[curve]);

    // TODO:
    // brushing
    // dragging genes around
    // rescaling axes

    /////////////////////////////////////////////////////////////////////////////
                          ///////    Set Up Chart    ///////
    /////////////////////////////////////////////////////////////////////////////

    var svg = d3.select("#graph-container").append("svg").attr("xmlns", "http://www.w3.org/2000/svg").attr("xmlns:xlink", "http://www.w3.org/1999/xlink");
    var g = svg.append("g").attr("transform", "translate("+margin.left+","+margin.top+")");
    svg.style("cursor", "move");

    var gene_set_name = "";
    var genes = [];
    var gene_wise = [];
    var sample_wise = [];
    var all_zero_genes = [];
    var x_scales = {};
    var genes_pc1 = [];
    var samples_pc1 = [];

    var x_scales;
    var title = g.selectAll(".title");
    var text = g.selectAll(".text");
    var axes = g.selectAll(".axis-x");
    var dots = g.selectAll(".dots");
    var line = g.selectAll(".line");
    var halo = g.selectAll(".halo");


    /////////////////////////////////////////////////////////////////////////////
                          ///////    Re-Draw Chart    ///////
    /////////////////////////////////////////////////////////////////////////////

    function restart({selected_gene_set_name=gene_set_name,
                      selected_genes=genes}={}) {

        gene_set_name = selected_gene_set_name;
        genes = selected_genes;
        matrix = _(samples_by_genes_matrix).mapObject((sample) => _.pick(sample, genes));

        // Gene-Wise
        gene_wise =
            _.zip(
                _(_(matrix).values()[0]).keys(),
                _.zip(..._(matrix).values().map((gene_counts) => _(gene_counts).values())).map((gene_counts) =>
                    _.object(_(matrix).keys(), gene_counts)
            )).map(([gene, counts]) => {
                counts_arr = Object.values(counts);
                all_zero = counts_arr.every(element => element === 0);
                mean = d3.mean(counts_arr);
                stddev = d3.deviation(counts_arr);
                mad = d3.mean(counts_arr.map((count) => Math.abs(mean - count)));
                return {
                'id'       : gene,
                'gene'     : gene,
                'all_zero' : all_zero,
                'min'      : d3.min(counts_arr),
                'max'      : d3.max(counts_arr),
                'mean'     : mean,
                'stddev'   : stddev,
                'mad'      : mad,
                'pc1'      : null,
                'counts'   : counts_arr,
                'samples'  : Object.values(_(counts).mapObject((count, sample) => { return {
                    'id'         : sample+"_"+gene,
                    'sample'     : sample,
                    'classes'    : classes[sample],
                    'gene'       : gene,
                    'count'      : count,
                    'zscore'     : (count - mean) / stddev,
                    'zscore_mad' : (count - mean) / mad,
                }})),
                'classes':
                    _(categories_to_values_to_members).mapObject((values_to_members) =>
                        _(values_to_members).mapObject((sample_ids, label) => {
                            class_counts = _(counts).pick(sample_ids);
                            class_mean = d3.mean(class_counts);
                            class_stddev = d3.deviation(class_counts);
                            class_mad = d3.mean(Object.values(class_counts).map((count) => Math.abs(class_mean - count)));
                            return {
                            'label'        : label,
                            'class_mean'   : class_mean,
                            'class_stddev' : class_stddev,
                            'class_mad'    : class_mad,
                            }
                        })
                    ),
                }
            });


        all_zero_genes = _(gene_wise.filter(gene => gene.all_zero)).pluck('gene');

        permutation_order = order();

        return {
            'gene_wise'   : gene_wise,
            'permutation_order': permutation_order,
        }

    }

    function order({values_='counts',sorting_='complete',show_all_zero_=show_all_zero}={}) {

        values = values_;
        sorting = sorting_;
        show_all_zero = show_all_zero_;

        // if (show_all_zero === false) { gene_wise = gene_wise.filter((gene) => gene.all_zero === false); }

        genes_pc1 = reorder.pca1d(reorder.transpose(gene_wise.map((gene) => gene.counts)));
        _(_.zip(gene_wise, genes_pc1)).each((gene, pc1) => gene.pc1 = pc1);

        if (gene_wise.length > 1) {
            // permutation_order = reorder.sort_order(genes_pc1);
            permutation_order = reorder.optimal_leaf_order()(gene_wise.map((gene) => gene.counts)).reverse();  // get dendogram out?
            gene_wise = reorder.stablepermute(gene_wise, permutation_order);

        } else { permutation_order = [0]; }

        // Sample-Wise
        sample_wise = _(Object.values(_(matrix).mapObject((genes, sample) => { return {
                                'id'     : sample,
                                'sample' : sample,
                                'classes': classes[sample],
                                'pc1'    : null,
                                'genes'  : reorder.stablepermute(Object.values(_(genes).mapObject((count, gene) => { return {
                                              'id'       : gene+"_"+sample,
                                              'gene'     : gene,
                                              'count'    : count,
                                              'all_zero' : _(all_zero_genes).contains(gene),
                                          }})),permutation_order)
                            }}))).flatten();

        // if (show_all_zero === false) { _(sample_wise).each((sample) => sample.genes = sample.genes.filter((gene) => gene.all_zero === false)); }

        samples_pc1 = reorder.pca1d(gene_wise.map((sample) => sample.counts));
        _(_.zip(sample_wise, samples_pc1)).each(([sample, pc1]) => { sample.pc1 = pc1; });
        samples_pc1_domain = d3.extent(samples_pc1);

        render();

        return permutation_order;

    }

    function render() {

        title = g.selectAll(".title");
        text = g.selectAll(".text");
        axes = g.selectAll(".axis-x");
        dots = g.selectAll(".dots");
        line = g.selectAll(".line");
        halo = g.selectAll(".halo");

        y = (index) => gene_spacing*index;

        // gene symbols
        text = text.data([]);
        text.exit().remove();
        text = text.data(gene_wise);
        text = text.enter()
                   .append("a")
                   .attr("class", "text")
                   .attr("xlink:href", (d) => "http://www.genecards.org/cgi-bin/carddisp.pl?gene="+d.gene)
                   .style("cursor", "pointer")
                       .append("text")
                       .attr("font-family", "sans-serif")
                       .text((d) => d.gene)
                       .style("text-anchor", "middle")
                       .attr("x", 0) // (x_pos(0) + x_neg(0)) / 2)
                       .attr("y", (d, i) => y(i))
                       .attr("dy", "1em");


        // scale for each gene
        x_scales = _.object(gene_wise.map((obj) =>
            [obj.gene, d3.scaleLinear().domain([0, d3.max([obj.mean*2, 10])]).range([50,h]).nice()]
        ));

        axes = axes.data([]);
        axes.exit().remove();  // can have a cute transition here
        axes = axes.data(Object.entries(x_scales))
                   .enter()
                   .append("g")
                   .attr("class", "axis axis-x")
                   .attr("transform", (d, i) => "translate(0," + y(i) + ")")
                   .each(function (d) { d3.select(this).call(d3.axisBottom(d[1])) });

        line = line.data([]);
        line.exit().remove();
        line = line.data(sample_wise)
                     .enter()
                     .append("path")
                     .attr("class", "line")
                     .style("fill", "none");

        dots = dots.data([]);
        dots.exit().remove();
        dots = dots.data(gene_wise)
                   .enter()
                   .append("g")
                   .attr("class", "dots")
                   .selectAll(".dot")
                       .data((d, i) => d.samples.map((sample) => Object.assign({'i':i}, sample)))
                       .enter()
                       .append("circle")
                       .attr("class", "dot")
                       .attr("cy", (d) => y(d.i) )
                       .attr("cx", (d) => x_scales[d.gene](d.count) );

        halo = halo.data([]);
        halo.exit().remove();
        halo = halo.data(sample_wise)
                     .enter()
                     .append("path")
                     .attr("class", "halo")
                     .style("fill", "none")
                     .style("stroke-width", 20)
                     .style("stroke-linecap", "round")
                     .style("stroke-linejoin", "round");

        title = title.data([]);
        title.exit().remove();
        title = title.data([gene_set_name])
                     .enter()
                     .append("text")
                     .attr("font-family", "sans-serif")
                     .attr("class", "title")
                     .text(gene_set_name)
                     .style("text-anchor", "middle")
                     .attr("x", 0)
                     .attr("y", 0)
                     .attr("dy", "-3em");

        style();

    }

    function style({show_points_=show_points,
                    point_radius_=point_radius,
                    point_coloring_system_=point_coloring_system,
                    points_color_by_=points_color_by,
                    default_point_color_=default_point_color,
                    show_lines_=show_lines,
                    show_halos_=show_halos,
                    line_coloring_system_=line_coloring_system,
                    lines_color_by_=lines_color_by,
                    default_line_color_=default_line_color,
                    negative_color_=negative_color,
                    middle_color_=middle_color,
                    positive_color_=positive_color,
                    curve_=curve,
                    halo_opacity_=halo_opacity}={}) {

        show_points = show_points_;
        point_radius = point_radius_;
        point_coloring_system = point_coloring_system_;
        points_color_by = points_color_by_;
        default_point_color = default_point_color_;
        show_lines = show_lines_;
        show_halos = show_halos_;
        line_coloring_system = line_coloring_system_;
        lines_color_by = lines_color_by_;
        default_line_color = default_line_color_;
        negative_color = negative_color_,
        middle_color = middle_color_,
        positive_color = positive_color_,
        curve = curve_;
        halo_opacity = halo_opacity_;

        if (line_coloring_system === 'z_stddev' || line_coloring_system === 'z_mad') {
            lines_color_by = gene_wise[0].gene;
            // bind click events?
        } else if (line_coloring_system === 'class') {
            lines_color_by = categories[0];
        } else if (line_coloring_system === 'pc1') {
            pc1_colors = d3.scaleLinear().domain(samples_pc1_domain).range([negative_color, positive_color]);
        }

        // dots are bound to gene_wise.samples
        svg.selectAll(".dot")
            .attr("visibility", show_points ? "visible" : "hidden")
            .attr("r", point_radius)
            .style("fill", (d) => point_colors[point_coloring_system](d));
            // .style("opacity", );


        // lines are bound to sample_wise
        line.attr("visibility", show_lines ? "visible" : "hidden")
            .attr("d", (d) => line_from_sample.curve(curves[curve])(d.genes))
            .style("stroke", (d) => line_colors[line_coloring_system](d));
             // .style("stroke-width", );
             // .style("stroke-opacity", );

        // halos are bound to sample_wise
        halo.attr("visibility", show_halos ? "visible" : "hidden")
            .attr("d", (d) => line_from_sample.curve(curves[curve])(d.genes))
            .style("stroke", (d) => line_colors[line_coloring_system](d))
            .style("stroke-opacity", halo_opacity);

    }


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
    }

}





