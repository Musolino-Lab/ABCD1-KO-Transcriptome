
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
    var line_bendiness = 0.5;

    var default_radius = 5;
    var gene_spacing = 30;

    // coloring
    // brushing

    /////////////////////////////////////////////////////////////////////////////
                          ///////    Set Up Chart    ///////
    /////////////////////////////////////////////////////////////////////////////

    var svg = d3.select("#graph-container").append("svg").attr("xmlns", "http://www.w3.org/2000/svg").attr("xmlns:xlink", "http://www.w3.org/1999/xlink");
    var g = svg.append("g").attr("transform", "translate("+margin.left+","+margin.top+")");
    svg.style("cursor", "move");

    var title = g.selectAll(".title");
    var axes = g.selectAll(".axis-x");
    var dots = g.selectAll(".dot");
    var lines = g.selectAll(".line");

    var gene_set_name = "";
    var genes = [];

    /////////////////////////////////////////////////////////////////////////////
                          ///////    Re-Draw Chart    ///////
    /////////////////////////////////////////////////////////////////////////////

    function restart({selected_gene_set_name=gene_set_name, selected_genes=genes}={}) {

        var gene_set_name = selected_gene_set_name;
        var genes = selected_genes;
        var matrix = _(samples_by_genes_matrix).mapObject((sample) => _.pick(sample, genes));

        var sample_wise = _(Object.values(_(matrix).mapObject((genes, sample) => { return {
                                'sample': sample,
                                'genes' : Object.values(_(genes).mapObject((count, gene) => { return {
                                                'gene': gene,
                                                'count': count,
                                            }}))
                            }}))).flatten();

        var gene_wise =
            _.zip(
                _(_(matrix).values()[0]).keys(),
                _.zip(..._(matrix).values().map((gene_counts) => _(gene_counts).values())).map((gene_counts) =>
                    _.object(_(matrix).keys(), gene_counts)
            )).map(([gene, counts]) => {
                counts_arr = Object.values(counts),
                mean = d3.mean(counts_arr);
                stddev = d3.deviation(counts_arr);
                mad = d3.mean(counts_arr.map((count) => Math.abs(mean - count)));
                return {
                'gene': gene,
                'mean': mean,
                'stddev': stddev,
                'mad': mad,
                'counts': counts_arr,
                'samples': Object.values(_(counts).mapObject((count, sample) => { return {
                    'sample': sample,
                    'gene': gene,
                    'count': count,
                    'zscore': (count - mean) / stddev,
                    'zscore_mad': (count - mean) / mad,
                }})),
                'classes': Object.values(_(classes).mapObject((sample_ids, label) => {
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

        console.log(gene_wise);


        var data = _(Object.values(_(matrix).mapObject((genes, sample) =>
                        Object.values(_(genes).mapObject((count, gene) => {return {
                             'sample': sample,
                             'gene': gene,
                             'count': count,
                        }}))))).flatten();
                    // .filter(not_null);

        var x_scales = _(_.object(
            _(_(matrix).values()[0]).keys(),
            _.zip(..._(matrix).values().map((gene_counts) => _(gene_counts).values())))
        ).mapObject((counts) => d3.scaleLinear().domain([0,d3.max(counts)*1.5]).range([0,h]).nice());

        // sort genes in some manner  TODO


        // create a scale object for all the genes
        var y = d3.scalePoint().domain(genes).range([0,gene_spacing*genes.length]);

        var y_axis = g.append("g")
                      .attr("class", "axis axis-y")
                      .call(d3.axisLeft(y));

        // draw a scale for each gene
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
        dots = dots.data(data)
        dots.enter()
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
                                  .curve(d3.curveCatmullRom.alpha(line_bendiness));

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

    }

}





