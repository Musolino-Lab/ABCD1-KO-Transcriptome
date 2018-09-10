

let range = n => [...Array(n).keys()];
let clamp = (min, max) => ((x) => Math.min(Math.max(x, min), max));
let transpose = (array) => array[0].map((col, i) => array.map(row => row[i]));
let flatten = (array) => [].concat.apply([], array);
let safeStr = (str) => str.split(' (')[0].replace(/\ /gi, '_');
let sum_counts_objects = (a, b) => _.object(_.uniq(Object.keys(a).concat(Object.keys(b))).map(key => [key, (a[key] || 0) + (b[key] || 0)]));

function Heatmap(samples_by_genes_matrix, gene_sets, classes, separate_by) {

    // samples_by_genes_matrix: {'sample_id': {'gene1': float, 'gene2': float}};   // since this is an object, sample IDs and gene names are guaranteed to be unique.

    // gene_sets: {'gene_set_name': ['gene1', 'gene2']};

    // classes: {'sample_id': {'category1': 'value', 'category2': value}}

    // categories: ['column_name_of_sample_class_labels_1', 'column_name_of_sample_class_labels_2']
    var categories = _.uniq(_.flatten(Object.values(classes).map(obj => Object.keys(obj))));

    // categories_to_members_to_values: {'category': {'sample_id1': 'value', 'sample_id2': value}}
    var categories_to_members_to_values = _.object(categories.map((category) =>
        [category, _.object(Object.entries(classes).map(([sample_id, categories_to_values]) => [sample_id, categories_to_values[category]]))]
    ));

    // categories_to_values_to_members: {'category': {'value1': ['sample_id1', 'sample_id2']}}
    var categories_to_values_to_members = _(categories_to_members_to_values).mapObject((samples_to_values) =>
        _(_(Object.entries(samples_to_values)).groupBy(([sample_id, value]) => value)).mapObject((arr_of_arr) => arr_of_arr.map(arr => arr[0]))
    );

    var hierarchy = {'id':'metadata', 'children':[]};
    Object.entries(classes).forEach(([sample_id, metadata]) => {
        pointer = hierarchy.children;
        prefix = 'metadata';
        categories.forEach((category, i) => {
            value = metadata[category];
            prefix += '-'+safeStr(value);
            existing_index_for_value = _(pointer).findIndex({'id':prefix});
            if (existing_index_for_value > -1) {
                if (i+1 === categories.length) { pointer[existing_index_for_value].children.push({'id':prefix+"-"+sample_id, 'name':sample_id}); }
                else { pointer = pointer[existing_index_for_value].children; }
            } else {
                if (i+1 === categories.length) { pointer.push({'id':prefix,'name':value,'category':category,'children':[{'id':prefix+"-"+sample_id, 'name':sample_id}]}); }
                else {
                    pointer.push({'id':prefix,'name':value,'category':category,'children':[]});
                    pointer = pointer[_(pointer).findIndex({'id':prefix})].children;
                }
            }
        })
    });

    var metadata = d3.hierarchy(hierarchy);
    var sample_to_sample_id = _.object(metadata.leaves().map(leaf => [leaf.data.name, leaf.data.id]));

    var internal_node_ordering = _(categories_to_values_to_members).mapObject(values_to_members => _.object(Object.keys(values_to_members).map((value, i) => [value, i])));
    metadata.each(node => { if (node.height !== 0 && node.depth !== 0) { node.data.order = internal_node_ordering[node.data.category][node.data.name]; } });

    var samples_to_bin = _(classes).mapObject(categories_to_values => Object.entries(_(categories_to_values).pick(separate_by)).sort().reduce((acc, [category, value]) => (acc ? acc+'-'+value : value), ''));
    var bin_to_samples = _(Object.keys(samples_to_bin)).groupBy(sample => samples_to_bin[sample]);

    var samples_to_small_bin = _(classes).mapObject(categories_to_values => Object.entries(categories_to_values).sort().reduce((acc, [category, value]) => (acc ? acc+'-'+value : value), ''));
    var small_bin_to_samples = _(Object.keys(samples_to_small_bin)).groupBy(sample => samples_to_small_bin[sample]);


    /////////////////////////////////////////////////////////////////////////////
                    ///////    Structure Variables    ///////
    /////////////////////////////////////////////////////////////////////////////

    var margin = {top: 100, right: 100, bottom: 0, left: 100};

    var w = window.innerWidth - (margin.left + margin.right);
    var h = window.innerHeight - (margin.top + margin.bottom);

    var values = 'zscore_stddev';
    var show_averages = false;
    var sorting = 'complete';
    var reordering = true;
    var minimum_nonzero = 0;

    value_accessors = {
        counts: (gene) => gene.counts,
        zscore_stddev: (gene) => gene.samples.map((sample) => sample.zscore_stddev),
        zscore_mad: (gene) => gene.samples.map((sample) => sample.zscore_mad),
        pc1: (gene) => gene.samples.map((sample) => sample.pc1),
    }
    value_accessor = value_accessors.counts;

    var ordered_sample_ids = [];
    var ordered_gene_ids = [];

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

    var category_colors = _.object(categories.map((category) => [category, d3.scaleOrdinal(d3.schemeCategory10)]))

    var colors = {};

    var spacing = 2;
    var rect_width = 16;
    var rect_height = 16;

    /////////////////////////////////////////////////////////////////////////////
                          ///////    Set Up Chart    ///////
    /////////////////////////////////////////////////////////////////////////////

    var svg = d3.select("#graph-container").append("svg").attr("xmlns", "http://www.w3.org/2000/svg").attr("xmlns:xlink", "http://www.w3.org/1999/xlink");
    var g = svg.append("g");
    svg.style("cursor", "move");

    var selected_gene_sets = {}
    var genes = [];
    var gene_wise = [];
    var gene_wise_indexer = {};
    var ordered_gene_wise = [];
    var sample_wise = [];
    var sample_wise_indexer = {};

    var legend_color = g.append("g").attr("class", "legend legendColor");


    /////////////////////////////////////////////////////////////////////////////
                          ///////    Methods    ///////
    /////////////////////////////////////////////////////////////////////////////

    function restart({selected_gene_sets_=selected_gene_sets}={}) {

        selected_gene_sets = selected_gene_sets_;
        // TODO deal with duplicate gene sets E.g. [Peroxisome, Peroxisome]
        keys = Object.keys(Object.values(samples_by_genes_matrix)[0]);  // genes included in matrix
        genes = d3.hierarchy({
                    'id': 'genes',
                    'children': selected_gene_sets.map((gs, i) => { return {
                        'id': safeStr(gs.gene_set_name),
                        'children': gs.genes.filter(gene => keys.includes(gene)).map(gene => {
                             return {'gene_set':gs.gene_set_name, 'gene':gene, 'id':(safeStr(gs.gene_set_name) || 'other')+'_'+gene}
                        })
                    }})
                });


        matrix = _(samples_by_genes_matrix).mapObject((sample) => _(sample).pick(genes.leaves().map(d => d.data.gene)));

        sample_wise = Object.entries(matrix).map(([sample, genes]) =>
            Object.entries(genes).map(([gene, count]) => { return {
                'id'        : sample+"_"+gene,
                'sample'    : sample,
                'sample_id' : sample_to_sample_id[sample],
                'gene'      : gene,
                'count'     : count,
            }})
         );

        gene_wise = transpose(sample_wise).map((by_gene) => {
            bin = _.object(Object.entries(bin_to_samples).map(([bin, samples]) => {
                l = by_gene.filter(gene => samples.includes(gene.sample));
                min    = d3.min(l, gene => gene.count);
                max    = d3.max(l, gene => gene.count);
                mean   = d3.mean(l, gene => gene.count);            // should we only be counting non-zeros?
                stddev = d3.deviation(l, gene => gene.count);
                mad    = d3.mean(l.map(gene => Math.abs(mean - gene.count)));
                return [bin, {'min':min,'max':max,'mean':mean,'stddev':stddev,'mad':mad}];
            }));
            num_nonzeros = by_gene.filter(d => d.count !== 0).length;
            return by_gene.map((gene) => Object.assign(gene, {
                'min'      : bin[samples_to_bin[gene.sample]].min,
                'max'      : bin[samples_to_bin[gene.sample]].max,
                'mean'     : bin[samples_to_bin[gene.sample]].mean,
                'stddev'   : bin[samples_to_bin[gene.sample]].stddev,
                'mad'      : bin[samples_to_bin[gene.sample]].mad,
                'num_nonzeros'  : num_nonzeros,
                'zscore_stddev' : bin[samples_to_bin[gene.sample]].stddev === 0 ? 0 : (gene.count - bin[samples_to_bin[gene.sample]].mean) / bin[samples_to_bin[gene.sample]].stddev,
                'zscore_mad'    : bin[samples_to_bin[gene.sample]].mad === 0 ? 0    : (gene.count - bin[samples_to_bin[gene.sample]].mean) / bin[samples_to_bin[gene.sample]].mad,
            }));
        });

        gene_wise_indexer =  _.object(gene_wise.map((gene, i) => [gene[0].gene, i]));
        sample_wise_indexer = _.object(sample_wise.map((sample, i) => [sample[0].sample, i]));

        order();
        render();

    }

    function order({values_=values,
                    sorting_=sorting,
                    minimum_nonzero_=minimum_nonzero}={}) {

        values = values_;
        value_accessor = value_accessors[values];
        sorting = sorting_;
        minimum_nonzero = minimum_nonzero_;

        // Filter by number non-zeros
        ordered_gene_wise = genes.leaves().map(leaf => gene_wise[gene_wise_indexer[leaf.data.gene]].map(sample => Object.assign({'gene_id':leaf.data.id}, sample)));

        ordered_gene_wise = ordered_gene_wise.filter((gene) => gene[0].num_nonzeros >= minimum_nonzero);

        if (ordered_gene_wise.length === 0) { return; }  // do something smart here.

        // Set Gene-Wise PC1
        // genes_pc1 = reorder.pca1d(reorder.transpose(ordered_gene_wise.map(value_accessor)));
        // _(_.zip(ordered_gene_wise, genes_pc1)).each((gene, pc1) => gene.pc1 = pc1);

        if (reordering && ordered_gene_wise.length > 1) {

            genes.each(node => {
                if (node.depth === 1 && node.name !== null) {  // for each gene set

                    set_gene_ids = node.children.map(gene => gene.data.id);
                    set_genes = ordered_gene_wise.filter(bygene => set_gene_ids.includes(bygene[0].gene_id));
                    set_gene_values = set_genes.map(bygene => bygene.map(gene => gene[values]));

                    if (sorting === 'complete') { permutation_order = reorder.optimal_leaf_order()(set_gene_values); } // get dendogram out?
                    else if (sorting === 'pc1') { permutation_order = reorder.sort_order(genes_pc1); }

                    set_gene_ids = reorder.stablepermute(set_gene_ids, permutation_order);
                    permutation_order = _.object(set_gene_ids.map((id, i) => [id, i]));
                    node.children.forEach((gene) => gene.data.order = permutation_order[gene.data.id]);

                }
            });

            Object.entries(small_bin_to_samples).forEach(([bin, samples]) => {

                samples = sample_wise.filter(by_sample => samples.includes(by_sample[0].sample));
                sample_ids = samples.map(by_sample => by_sample[0].sample_id);
                sample_values = samples.map(by_sample => by_sample.map(sample => sample[values]));

                permutation_order = reorder.optimal_leaf_order()(sample_values);
                sample_ids = reorder.stablepermute(sample_ids, permutation_order);
                permutation_order = _.object(sample_ids.map((id, i) => [id, i]));
                metadata.leaves().forEach((sample) => { if (sample.data.id in permutation_order) {sample.data.order = permutation_order[sample.data.id]} });

            });

        }

        genes = genes.count().sort(function(a, b) { return b.height - a.height || a.data.order - b.data.order; });
        ordered_gene_ids = genes.leaves().map(leaf => leaf.data.id);

        metadata = metadata.count().sort(function(a, b) { return b.height - a.height || a.data.order - b.data.order; });
        ordered_sample_ids = metadata.leaves().map(leaf => leaf.data.id);

    }

    function render({}={}) {

        margin = {'2': 10};

        metadata_across = rect_width*(categories.length+2);
        metadata_topdown = rect_height*ordered_sample_ids.length+2;
        d3.partition().size([metadata_topdown, metadata_across]).padding(2).round(true)(metadata);

        genes_across = (rect_height*ordered_gene_ids.length)+2;
        genes_topdown = 200;
        d3.partition().size([genes_across, genes_topdown]).padding(2).round(true)(genes);

        counter = 0; current_depth = 0;
        metadata.each((node) => {
            node.num_below = _(node.descendants().map(d => d.height)).countBy();

            if (node.depth !== current_depth) { current_depth = node.depth; counter = 0; }
            node.sibling_index = counter; counter += 1;

            node.offset = 0;
            if (node.parent) {
                node.num_left = node.parent.children.filter(sibling => sibling.sibling_index !== undefined && sibling.sibling_index < node.sibling_index)
                                                    .reduce((acc, sibling) => sum_counts_objects(acc, sibling.num_below), {});

                node.offset = d3.sum(Object.entries(node.num_left).map(([level, num]) => (margin[level] || 0)*(num))) + node.parent.offset;
                node.x0 += node.offset;
                node.x1 += node.offset + d3.sum(Object.entries(node.num_below).map(([level, num]) => (margin[level] || 0)*(num-1)));
            }

        });


        x = _.object(metadata.leaves().map(leaf => [leaf.data.id, leaf.x0]))
        y = _.object(genes.leaves().map(leaf => [leaf.data.id, leaf.x0]))

        category_y = (category) => -20*categories.indexOf(category);

        rect = g.selectAll(".rect").data(flatten(ordered_gene_wise), d => d.id);
        gene = g.selectAll(".gene").data(genes.descendants(), d => d.data.id);
        meta = g.selectAll(".meta").data(metadata.descendants(), d => d.data.id);

        // phase 1
            // rectangles which are exiting fade out
            // gene names which are exiting fade out
            // gene groups which are exiting fade out
        t_last = d3.transition().duration(200);
        if (rect.exit().size() > 0) {
            rect.exit().transition(t_last).style("opacity", 0).remove();
            gene.exit().transition(t_last).style("opacity", 0).remove();
            meta.exit().transition(t_last).style("opacity", 0).remove();
            t_last = t_last.transition().duration(500);
        }

        // phase 2
            // re-arrange ROWS (rectangles, gene symbols, gene groups)
        rect.transition(t_last).attr('y', d => y[d.gene_id])
        gene.transition(t_last).attr('y', d => y[d.data.id])
        meta.filter(function(d) { return d.children !== undefined && d.depth > 0; }).transition(t_last).attr('y', d => d.y0 - (categories.length+1)*rect_height - 8)
        meta.filter(d => d.height === 0).transition(t_last).attr('transform', function (d) {
            current_x = d3.select(this).attr('transform').split("translate(")[1].split(",")[0];
            return "translate("+current_x+","+ordered_gene_ids.length*rect_height+")rotate(60)" });
        t_last = t_last.transition().duration(500);

        // phase 3
            // re-arrange COLUMNS (rectangles, sample names, meta)
        rect.transition(t_last).attr('x', d => x[d.sample_id]);
        // meta.filter(d => d.height === 0).transition(t_last); change transform here
        meta.filter(function(d) { return d.children !== undefined && d.depth > 0; }).transition(t_last).attr('x', d => d.x0);
        t_last = t_last.transition().duration(500);

        // phase 4
            // rectangles which are entering get appended
            // gene names which are entering get appended
            // gene groups which are entering get appended
        rect.enter()
            .append('rect')
            .attr('class', 'rect')
            .attr('id', d => d.id)
            .attr('x', d => x[d.sample_id])
            .attr('y', d => y[d.gene_id])
            .attr('width', rect_width-2)
            .attr('height', rect_height-2)
            // .attr('fill', d => colors[values](d[values]))
            .style("opacity", 0).transition(t_last).style("opacity", 1);


        gene.enter()
            .filter(d => d.height === 0)  // leaves
            .append('text')
            .attr('class', 'gene')
            .attr('id', d => d.data.id)
            .attr('x', -10)
            .attr('y', d => d.x0)
            .text(d => d.data.gene)
            .attr("font-family", "sans-serif")
            .style("font-weight", 300)
            .style("cursor", "pointer")
            .style("text-anchor", "end")
            .attr("dy", "0.8em")
            .on("click", (d) => GeneCards(d.data.gene))
            .call(d3.drag().on("start", drag_gene_start).on("drag", drag_gene).on("end", drag_gene_end))
            .style("opacity", 0).transition(t_last).style("opacity", 1);

        gene.enter()
            .filter(function(d) { return d.depth === 1; })  // gene sets
            .append('rect')
            .attr('class', 'gene')
            .attr('id', d => d.data.id)
            .attr('x', d => -150)
            .attr('y', d => d.x0)
            .attr("width", function(d) { return d.y1 - d.y0; })
            .attr("height", function(d) { return d.x1 - d.x0; })
            .style('fill-opacity', 0)
            .style('stroke', 'black')
            .call(d3.drag().on("start", drag_gs_start).on("drag", drag_gs).on("end", drag_gs_end))
            .style("opacity", 0).transition(t_last).style("opacity", 1);


        meta.enter()
            .filter(d => d.height === 0)  // leaves
            .append('text')
            .attr('class', 'meta')
            .attr('id', d => d.data.id)
            .attr('transform', d => "translate("+(d.x0+rect_width/2)+","+ordered_gene_ids.length*rect_height+")rotate(60)")
            .text(d => d.data.name)
            .attr("font-family", "sans-serif")
            .style("font-weight", 300)
            .style("cursor", "pointer")
            .style("text-anchor", "start")
            .attr("dy", "0.8em")
            .call(d3.drag().on("start", drag_sample_start).on("drag", drag_sample).on("end", drag_sample_end))
            .style("opacity", 0).transition(t_last).style("opacity", 1);

        meta.enter()
            .filter(function(d) { return d.children !== undefined && d.depth > 0; })  // internal nodes
            .append("g")
            .attr("class", "meta")
            .attr("id", d => d.data.id)
            .attr('transform', d => 'translate('+d.x0+','+(d.y0 - (categories.length+1)*rect_height - 8)+')')
            // .attr('transform', d => 'translate('+d.x0+','+category_y(d.data.category)+')')
            .call(d3.drag().on("start", drag_meta_start).on("drag", drag_meta).on("end", drag_meta_end))
                .append('rect')
                .attr('class', 'category_box')
                .attr("width", function(d) { return d.x1 - d.x0; })
                .attr("height", function(d) { return d.y1 - d.y0; })
                .style("fill", d => category_colors[d.data.category](d.data.name))
            .select(function() { return this.parentNode; })
                .append('text')
                .attr('class', 'category_label')
                .text(d => d.data.name)
                .attr("font-family", "sans-serif")
                .style("text-anchor", "start")
                .style("font-size", 10)
                .attr("dy", "1.2em")
                .attr("dx", "0.2em")
            .select(function() { return this.parentNode; })
            .style("opacity", 0).transition(t_last).style("opacity", 1);

          // cell.append("clipPath")
          //     .attr("id", function(d) { return "clip-" + d.id; })
          //   .append("use")
          //     .attr("xlink:href", function(d) { return "#rect-" + d.id + ""; });

          // cell.append("text")
          //     .attr("clip-path", function(d) { return "url(#clip-" + d.id + ")"; })



    style();

    }

    function style({negative_color_=negative_color,
                    middle_color_=middle_color,
                    positive_color_=positive_color,
                    show_legends_=show_legends}={}) {

        negative_color = negative_color_,
        middle_color = middle_color_,
        positive_color = positive_color_,
        show_legends=show_legends_;

        all_values = flatten(gene_wise);

        colors = {
            zscore_stddev:  d3.scaleLinear().domain([d3.min(all_values, d => d.zscore_stddev), 0, d3.max(all_values, d => d.zscore_stddev)]).range([negative_color, middle_color, positive_color]),
            zscore_mad:     d3.scaleLinear().domain([d3.min(all_values, d => d.zscore_mad),    0, d3.max(all_values, d => d.zscore_mad)]   ).range([negative_color, middle_color, positive_color]),
            pc1:            d3.scaleLinear().domain([d3.min(all_values, d => d.pc1),           0, d3.max(all_values, d => d.pc1)]          ).range([negative_color, middle_color, positive_color]),
        }

        g.selectAll(".rect")
            .style("fill", (d) => colors[values](d[values]));



        if (show_legends) {

        } else {

        }

    }


    /////////////////////////////////////////////////////////////////////////////
                          ///////    Drag Axes    ///////
    /////////////////////////////////////////////////////////////////////////////

    function drag_gene_start(d) {


        // console.log(ordered_gene_wise);

    }

    function drag_gene(d) {

        dragged_index = _(ordered_gene_wise).findIndex((gene) => gene.gene === d);
        // console.log(dragged_index);

        // g.select("#"+d).attr("transform", function(d, i) { return "translate(0," +  expr(i) + ")" });
        // g.selectAll(".rect").attr("transform", function(d, i) { return "translate(0," + (expr(i) - max_point_radius) + ")" });


        // let expr = (current_index) => {
        //     if (current_index < dragged_index) {
        //         if (current_index < ((d3.event.y-s()/2) / s())) {
        //                 return y(current_index);
        //             } else {
        //                 return y(current_index) + s();
        //             }
        //     } else {  // current_index >= dragged_index
        //         if (current_index < ((d3.event.y+s()/2) / s())) {
        //                 return y(current_index) - s();
        //             } else {
        //                 return y(current_index);
        //             }
        //     }
        // }


        // d3.select(this).attr("transform", "translate(0," + d3.event.y + ")");
        // d3.select(".dots#"+d[0]).attr("transform", (d, i) => "translate(0," + (d3.event.y - max_point_radius) + ")");

    }

    function drag_gene_end(d) {

        // dragged_index = _(ordered_gene_wise).findIndex((gene) => gene.id === d[0]);
        // old_index = dragged_index;
        // new_index = clamp(0, ordered_gene_wise.length)(Math.round(d3.event.y / s()));

        // ordered_gene_wise.splice(new_index, 0, ordered_gene_wise.splice(old_index, 1)[0]);
        // sample_wise.forEach((sample) => sample.genes.splice(new_index, 0, sample.genes.splice(old_index, 1)[0]));

        // render();

    }


    function drag_sample_start(d) {
        // console.log(sample_wise);

    }
    function drag_sample(d) {}
    function drag_sample_end(d) {}


    function drag_meta_start(d) {
        console.log(d);

    }
    function drag_meta(d) {}
    function drag_meta_end(d) {}

    function drag_gs_start(d) {
        console.log(d);

    }
    function drag_gs(d) {}
    function drag_gs_end(d) {}

    /////////////////////////////////////////////////////////////////////////////
                          ///////    Brush Axes    ///////
    /////////////////////////////////////////////////////////////////////////////


    function setFocus(d) {
    }

    function removeFocus(d) {
    }

    function GeneCards(d) { window.open("http://www.genecards.org/cgi-bin/carddisp.pl?gene="+d,'_blank') }


    /////////////////////////////////////////////////////////////////////////////
                          ///////   Zoom & Resize    ///////
    /////////////////////////////////////////////////////////////////////////////

    svg.call(d3.zoom().on("zoom", zoomed)).on("wheel.zoom", wheeled);

    transform = d3.zoomTransform(g);
    transform.x += margin.left;
    transform.y += margin.top;
    g.attr("transform", transform);

    function zoomed() {
        current_transform = d3.zoomTransform(g);
        current_transform.x += d3.event.sourceEvent.movementX;
        current_transform.y += d3.event.sourceEvent.movementY;
        g.attr("transform", current_transform);
    }

    function wheeled() {
        current_transform = d3.zoomTransform(g);
        if (d3.event.ctrlKey) {
            current_transform.k = clamp(0.1, 5)(current_transform.k - d3.event.deltaY * 0.01);
        } else {
            current_transform.y = clamp(-(ordered_gene_ids.length*rect_height-100), h)(current_transform.y - d3.event.deltaY);
        }
        g.attr("transform", current_transform);
    }

    function resize() {
        svg.attr("width", $("#graph-container").innerWidth()).attr("height", $("#graph-container").innerHeight());
        w = $("#graph-container").innerWidth() - (margin.left + margin.right);
        h = $("#graph-container").innerHeight() - (margin.top + margin.bottom);
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

        get_sorted_gene_list: () => _(sample_wise[0]).pluck('gene'),

        set_reordering: (reordering_) => { reordering = reordering_; if (reordering) { order(); render(); } },
    }

}





