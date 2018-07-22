
function Braid(samples_by_genes_matrix, gene_sets, classes) {



    /////////////////////////////////////////////////////////////////////////////
                          ///////    Variables    ///////
    /////////////////////////////////////////////////////////////////////////////

    var margin = {top: 100, right: 100, bottom: 0, left: 100};

    var w = window.innerWidth - (margin.left + margin.right);
    var h = window.innerHeight - (margin.top + margin.bottom);

    var absolute_axis = true;
    var relative_axis = false;
    var margin_between_axes = 10;

    var show_lines_connecting_samples = true;
    var line_bendiness = 0;

    var default_radius = 3;
    var gene_spacing = 50;

    // coloring  // dot with the first PC?
    // brushing

    var classes_to_members = _(Object.entries(classes)).groupBy(sample_and_class => sample_and_class[1]);

    // var sample_colors = ;
    // var class_colors = ;

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

    var title = g.selectAll(".title");
    var axes = g.selectAll(".axis-x");
    var dots = g.selectAll(".dots");
    var lines = g.selectAll(".line");

    /////////////////////////////////////////////////////////////////////////////
                          ///////    Re-Draw Chart    ///////
    /////////////////////////////////////////////////////////////////////////////

    function restart({selected_gene_set_name=gene_set_name, selected_genes=genes}={}) {

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
                'id': gene,
                'gene': gene,
                'all_zero': all_zero,
                'min': d3.min(counts_arr),
                'max': d3.max(counts_arr),
                'mean': mean,
                'stddev': stddev,
                'mad': mad,
                'counts': counts_arr,
                'samples': Object.values(_(counts).mapObject((count, sample) => { return {
                    'id' : sample+"_"+gene,
                    'sample': sample,
                    'class': classes[sample],
                    'gene': gene,
                    'count': count,
                    'zscore': (count - mean) / stddev,
                    'zscore_mad': (count - mean) / mad,
                }})),
                'classes': Object.values(_(classes_to_members).mapObject((sample_ids, label) => {
                    class_counts = _(counts).pick(sample_ids);
                    class_mean = d3.mean(class_counts);
                    class_stddev = d3.deviation(class_counts);
                    class_mad = d3.mean(Object.values(class_counts).map((count) => Math.abs(class_mean - count)));
                    return {
                    'label': label,
                    'class_mean': class_mean,
                    'class_stddev': class_stddev,
                    'class_mad': class_mad,
                    }
                })),
            }});

        all_zero_genes = _(gene_wise.filter(gene => gene.all_zero)).pluck('gene');

        genes_pc1 = reorder.pca1d(reorder.transpose(gene_wise.map((gene) => gene.counts)));
        _(_.zip(gene_wise, genes_pc1)).each((gene, pc1) => gene.pc1 = pc1);

        // Ordering
        // perm = reorder.sort_order(genes_pc1);
        perm = reorder.optimal_leaf_order()(gene_wise.map((gene) => gene.counts)).reverse();  // get dendogram out?
        gene_wise = reorder.stablepermute(gene_wise, perm);

        // Sample-Wise
        sample_wise = _(Object.values(_(matrix).mapObject((genes, sample) => { return {
                                'id': sample,
                                'sample': sample,
                                'class' : classes[sample],
                                'genes' : reorder.stablepermute(Object.values(_(genes).mapObject((count, gene) => { return {
                                              'id': gene+"_"+sample,
                                              'gene': gene,
                                              'count': count,
                                              'all_zero': _(all_zero_genes).contains(gene),
                                          }})),perm)
                            }}))).flatten();

        samples_pc1 = reorder.pca1d(gene_wise.map((sample) => sample.counts));
        _(_.zip(sample_wise, samples_pc1)).each((sample, pc1) => sample.pc1 = pc1);

        return render();
    }

    function render() {

        // scale for all genes
        var y = d3.scalePoint().domain(gene_wise.map(o => o.gene)).range([0,gene_spacing*genes.length]);
        var y_axis = g.append("g").attr("class", "axis axis-y").call(d3.axisLeft(y));

        // scale for each gene
        x_scales = _.object(gene_wise.map((obj) =>
            [obj.gene, d3.scaleLinear().domain([0, obj.mean*2]).range([0,h]).nice()]
        ));

        axes = axes.data([]);
        axes.exit().remove();  // can have a cute transition here
        axes = axes.data(Object.entries(x_scales));
        axes.enter()
            .append("g")
            .attr("class", "axis axis-x")
            .attr("transform", (d) => "translate(0," + y(d[0]) + ")")
            .each(function (d) { d3.select(this).call(d3.axisBottom(d[1])) });

        // dots
        dots = dots.data([]);
        dots.exit().remove();
        dots = dots.data(gene_wise)
        dots.enter()
            .append("g")
            .attr("class", "dots")
            .selectAll(".dot")
                .data((d) => d.samples)
                .enter()
                .append("circle")
                .attr("class", "dot")
                .attr("r", default_radius)
                .attr("cy", (d) => y(d.gene) )
                .attr("cx", (d) => x_scales[d.gene](d.count) )
                .style("fill", "red");

        // lines
        var line_from_sample = d3.line()
                                  .y((d) => y(d.gene))
                                  .x((d) => x_scales[d.gene](d.count))
                                  .curve(d3.curveLinear);
                                  // .curve(d3.curveCatmullRom.alpha(line_bendiness));
                                  // .curve(d3.curveNatural);
                                  // .curve(d3.curveMonotoneY);
                                  // .curve(d3.curveCardinal.tension(1-line_bendiness));

        lines = lines.data([]);
        lines.exit().remove()
        if (show_lines_connecting_samples) {
        lines = lines.data(sample_wise)
        lines.enter()
             .append("path")
             .attr("d", (d) => line_from_sample(d.genes))
             .attr("stroke", "blue")
             .style("fill", "none");
        }

        // title
        title = title.data([]);
        title.exit().remove();
        title = title.data([gene_set_name]);
        title = title.enter()
                     .append("text")
                     .attr("font-family", "sans-serif")
                     .attr("class", "title")
                     .text(gene_set_name)
                     .style("text-anchor", "middle")
                     .attr("x", 0)
                     .attr("y", 0)
                     .attr("dy", "-3em");

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
        'restart': restart,
        'render' : render,
        'gene_wise': gene_wise,
        'sample_wise': sample_wise,

    }

}





