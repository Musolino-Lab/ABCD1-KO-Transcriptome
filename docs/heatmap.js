

let range = n => [...Array(n).keys()];
let clamp = (min, max) => ((x) => Math.min(Math.max(x, min), max));
let transpose = (array) => array[0].map((col, i) => array.map(row => row[i]));
let flatten = (array) => [].concat.apply([], array);
let safeStr = (str) => str.split(' (')[0].replace(/\ /gi, '_');
let sum_counts_objects = (a, b) => _.object(_.uniq(Object.keys(a).concat(Object.keys(b))).map(key => [key, (a[key] || 0) + (b[key] || 0)]));
let pointing_right = (d) => ''+((d.x1 - d.x0) * 1 + (d.y1 - d.y0) * 1)+','+(d.x1 - d.x0);  // https://stackoverflow.com/questions/8976791/how-to-set-a-stroke-width1-on-only-certain-sides-of-svg-shapes
let pointing_left = (d) => '0,'+((d.x1 - d.x0) * 1)+','+((d.x1 - d.x0) * 1 + (d.y1 - d.y0) * 2);
Array.prototype.last = function() { return this[this.length - 1]; };
Array.prototype.move = function(from, to) { this.splice(to, 0, this.splice(from, 1)[0]); };


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


    var colors20 = ['#3366cc', '#dc3912', '#ff9900', '#109618', '#990099', '#0099c6', '#dd4477',
                    '#66aa00', '#b82e2e', '#316395', '#994499', '#22aa99', '#aaaa11', '#6633cc',
                    '#e67300', '#8b0707', '#651067', '#329262', '#5574a6', '#3b3eac'];

    var category_colors = _.object(categories.map((category) => [category, d3.scaleOrdinal(d3.schemeCategory10)]))

    var colors = {};

    var spacing = 2;
    var rect_width = 16;
    var rect_height = 16;
    var margins = {'x': {'2': 10}, 'y': {'1': 10}};

    var y_axis_leaves_position = 'before';
    var y_axis_nodes_position = 'before';
    var x_axis_leaves_position = 'after';
    var x_axis_nodes_position = 'after';
    var y_axis_leaves_x, y_axis_nodes_x, category_y, x_axis_leaves_y;
    var x_axis_leaves_rotation = (x_axis_leaves_position === 'before') ? -60 : 60;


    /////////////////////////////////////////////////////////////////////////////
                          ///////    Set Up Chart    ///////
    /////////////////////////////////////////////////////////////////////////////

    var svg = d3.select('#graph-container').append('svg').attr('xmlns', 'http://www.w3.org/2000/svg').attr('xmlns:xlink', 'http://www.w3.org/1999/xlink');
    var g = svg.append('g');
    svg.style('cursor', 'move');

    var selected_gene_sets = {}
    var genes = [];
    var gene_wise = [];
    var gene_wise_indexer = {};
    var ordered_gene_wise = [];
    var sample_wise = [];
    var sample_wise_indexer = {};
    var metadata = {};
    var sample_to_sample_id = {};
    var x, y, category_y;


    var legend_color = g.append('g').attr('class', 'legend legendColor');


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
                        'name': gs.gene_set_name,
                        'children': gs.genes.filter(gene => keys.includes(gene)).map(gene => {
                             return {'gene_set':gs.gene_set_name, 'name':gene, 'id':(safeStr(gs.gene_set_name) || 'other')+'_'+gene}
                        })
                    }})
                });


        matrix = _(samples_by_genes_matrix).mapObject((sample) => _(sample).pick(genes.leaves().map(d => d.data.name)));

        sample_wise = Object.entries(matrix).map(([sample, genes]) =>
            Object.entries(genes).map(([gene, count]) => { return {
                'id'        : sample+'_'+gene,
                'sample'    : sample,
                'sample_id' : null,  // assigned later
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

        // METADATA

        // to order the levels: change the order of categories

        hierarchy = {'id':'metadata', 'children':[]};
        Object.entries(classes).forEach(([sample_id, metadata]) => {
            pointer = hierarchy.children;
            prefix = 'metadata';
            categories.forEach((category, i) => {
                value = metadata[category];
                prefix += '-'+safeStr(value);
                existing_index_for_value = _(pointer).findIndex({'id':prefix});
                if (existing_index_for_value > -1) {
                    if (i+1 === categories.length) { pointer[existing_index_for_value].children.push({'id':prefix+'-'+sample_id, 'name':sample_id}); }
                    else { pointer = pointer[existing_index_for_value].children; }
                } else {
                    if (i+1 === categories.length) { pointer.push({'id':prefix,'name':value,'category':category,'children':[{'id':prefix+'-'+sample_id, 'name':sample_id}]}); }
                    else {
                        pointer.push({'id':prefix,'name':value,'category':category,'children':[]});
                        pointer = pointer[_(pointer).findIndex({'id':prefix})].children;
                    }
                }
            })
        });

        metadata = d3.hierarchy(hierarchy);
        sample_to_sample_id = _.object(metadata.leaves().map(leaf => [leaf.data.name, leaf.data.id]));
        sample_wise.forEach(by_sample => by_sample.forEach(sample => sample.sample_id = sample_to_sample_id[sample.sample]));

        internal_node_ordering = _(categories_to_values_to_members).mapObject(values_to_members => _.object(Object.keys(values_to_members).map((value, i) => [value, i])));
        metadata.each(node => { if (node.height !== 0 && node.depth !== 0) { node.data.order = internal_node_ordering[node.data.category][node.data.name]; } });

        // GENES

        // Filter by number non-zeros
        ordered_gene_wise = genes.leaves().map(leaf => gene_wise[gene_wise_indexer[leaf.data.name]].map(sample => Object.assign({'gene_id':leaf.data.id}, sample)));

        ordered_gene_wise = ordered_gene_wise.filter((gene) => gene[0].num_nonzeros >= minimum_nonzero);

        if (ordered_gene_wise.length === 0) { return; }  // do something smart here.

        // Set Gene-Wise PC1
        // genes_pc1 = reorder.pca1d(reorder.transpose(ordered_gene_wise.map(value_accessor)));
        // _(_.zip(ordered_gene_wise, genes_pc1)).each((gene, pc1) => gene.pc1 = pc1);

        if (reordering && ordered_gene_wise.length > 1) {

            counter = 0;
            genes.each(node => {
                if (node.height === 1 && node.data.name !== null) {  // for each gene set
                    node.data.order = counter; counter += 1;

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

            // SAMPLES

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

        genes.count().sort(function(a, b) { return b.height - a.height || a.data.order - b.data.order; });
        ordered_gene_ids = genes.leaves().map(leaf => leaf.data.id);

        metadata.count().sort(function(a, b) { return b.height - a.height || a.data.order - b.data.order; });
        ordered_sample_ids = metadata.leaves().map(leaf => leaf.data.id);

    }

    function offset(partition, margin) {

        counter = 0; current_depth = 0;
        partition.each((node) => {
            node.num_below = _(node.descendants().map(d => d.height)).countBy();

            if (node.depth !== current_depth) { current_depth = node.depth; counter = 0; }
            node.sibling_index = counter; counter += 1;
        });

        partition.each((node) => {
            if (node.parent) {
                node.num_left = node.parent.children.filter(sibling => sibling.sibling_index !== undefined && sibling.sibling_index < node.sibling_index)
                                                    .reduce((acc, sibling) => sum_counts_objects(acc, sibling.num_below), {});

                node.offset = d3.sum(Object.entries(node.num_left).map(([level, num]) => (margin[level] || 0)*(num))) + node.parent.offset;
                node.x0 += node.offset;
                node.x1 += node.offset + d3.sum(Object.entries(node.num_below).map(([level, num]) => (margin[level] || 0)*(num-1)));
            } else {
                node.offset = 0;
                node.num_left = 0;
                node.x1 += d3.sum(Object.entries(node.num_below).map(([level, num]) => (margin[level] || 0)*(num-1)));
            }
        });

    }

    function render({}={}) {

        metadata_across = rect_width*ordered_sample_ids.length+2;
        metadata_topdown = rect_height*(categories.length+2);
        d3.partition().size([metadata_across, metadata_topdown]).padding(2).round(true)(metadata);
        offset(metadata, margins['x']);

        genes_across = (rect_height*ordered_gene_ids.length)+2;
        genes_topdown = rect_width*genes.height*2;
        d3.partition().size([genes_across, genes_topdown]).padding(2).round(true)(genes);
        offset(genes, margins['y']);

        x = _.object(metadata.leaves().map(leaf => [leaf.data.id, leaf.x0]))
        y = _.object(genes.leaves().map(leaf => [leaf.data.id, leaf.x0]))

        y_axis_leaves_x = (y_axis_leaves_position === 'before') ? -10 : metadata.x1 + 10;
        y_axis_nodes_x = (y_axis_nodes_position === 'before') ? ((y_axis_leaves_position === 'before') ? -100 : -rect_width-10) : metadata.x1 + ((y_axis_leaves_position === 'before') ? rect_width+10 : 100-rect_width);

        category_y = _.object((x_axis_nodes_position === 'before') ? categories.map((c, i) => [c, -((categories.length-i) * rect_height)-10]) : categories.map((c, i) => [c, genes.x1 + ((categories.length-i) * rect_height)]));
        x_axis_leaves_y = (x_axis_leaves_position === 'before') ? ((x_axis_nodes_position === 'before') ? -(categories.length) * rect_height : 0) : genes.x1 + ((x_axis_nodes_position === 'before') ? 10 : categories.length * rect_height + rect_height+6);
        x_axis_leaves_rotation = (x_axis_leaves_position === 'before') ? -60 : 60;

        rect = g.selectAll('.rect').data(flatten(ordered_gene_wise), d => d.id);
        gene = g.selectAll('.gene').data(genes.descendants(), d => d.data.id);
        meta = g.selectAll('.meta').data(metadata.descendants(), d => d.data.id);
        catg = g.selectAll('.catg').data(categories, d => d);

        // phase 1
            // rectangles which are exiting fade out
            // gene names which are exiting fade out
            // gene groups which are exiting fade out
        t_last = d3.transition().duration(200);
        if (rect.exit().size() > 0 || gene.exit().size() > 0 || meta.exit().size() > 0) {
            rect.exit().transition(t_last).style('opacity', 0).remove();
            gene.exit().transition(t_last).style('opacity', 0).remove();
            meta.exit().transition(t_last).style('opacity', 0).remove();
            t_last = t_last.transition().duration(500);
        }

        // phase 2
            // re-arrange ROWS
        rect.transition(t_last).attr('y', d => y[d.gene_id])
        gene.filter(node => node.height === 0).transition(t_last).attr('y', d => y[d.data.id])
        gene.filter(gene => gene.height === 1).transition(t_last).attr('transform', d => 'translate('+y_axis_nodes_x+','+d.x1+')rotate(-90)')
        t_last = t_last.transition().duration(500);

        // phase 3
            // re-arrange COLUMNS
        rect.transition(t_last).attr('x', d => x[d.sample_id]);
        meta.filter(node => node.height === 0).transition(t_last).attr('x', d => d.x0)
                                                                 .attr('y', x_axis_leaves_y)
                                                                 .attr('transform', d => 'rotate('+x_axis_leaves_rotation+','+d.x0+','+x_axis_leaves_y+')');
        meta.filter(node => node.depth > 0 && node.height > 0).transition(t_last).attr('transform', d => 'translate('+d.x0+','+category_y[d.data.category]+')');
        meta.filter(node => node.depth > 0 && node.height > 0).selectAll('.category_box').transition(t_last).attr('width', d => d.x1 - d.x0);
        catg.transition(t_last).attr('x', metadata.x1+10).attr('y', d => category_y[d]);
        t_last = t_last.transition().duration(500);


        var drag_gene     = (d) => drag_node(     d, genes, 'y', 'gene', '.gene');
        var drag_gene_end = (d) => drag_node_end( d, genes, 'y', ordered_gene_ids, ordered_gene_wise, sample_wise);

        var drag_meta     = (d) => drag_node(     d, metadata, 'x', 'sample', '.meta');
        var drag_meta_end = (d) => drag_node_end( d, metadata, 'x', ordered_sample_ids, sample_wise, ordered_gene_wise);

        // phase 4
            // rectangles which are entering get appended
            // gene names / groups which are entering get appended
            // sample names / groups which are entering get appended
        meta.enter()
            .filter(d => d.height === 0)  // leaves
            .append('text')
            .attr('class', 'meta')
            .attr('id', d => d.data.id)
            .attr('x', d => d.x0)
            .attr('y', x_axis_leaves_y)
            .attr('dx', (x_axis_leaves_position === 'before' ? rect_width-2 : 0))
            .attr('transform', d => 'rotate('+x_axis_leaves_rotation+','+d.x0+','+x_axis_leaves_y+')')
            .text(d => d.data.name)
            .attr('font-family', 'sans-serif')
            .style('font-weight', 300)
            .style('cursor', 'pointer')
            .style('text-anchor', 'start')
            .call(d3.drag().on('drag', drag_meta).on('end', drag_meta_end))
            .style('opacity', 0).transition(t_last).style('opacity', 1);

        meta.enter()
            .filter(d => d.children !== undefined && d.depth > 0)  // internal nodes
            .append('g')
            .attr('class', 'meta')
            .attr('id', d => d.data.id)
            .attr('transform', d => 'translate('+d.x0+','+category_y[d.data.category]+')')
            .call(d3.drag().on('drag', drag_meta).on('end', drag_meta_end))
                .append('rect')
                .attr('class', 'category_box')
                .attr('id', d => 'rect-' + d.data.id)
                .attr('width', d => d.x1 - d.x0)
                .attr('height', d => d.y1 - d.y0)
                .style('fill', d => category_colors[d.data.category](d.data.name))
            .select(function() { return this.parentNode; })
                .append('text')
                .attr('class', 'category_label')
                .attr('clip-path', d => 'url(#clip-' + d.data.id + ')')
                .text(d => d.data.name)
                .attr('font-family', 'sans-serif')
                .style('text-anchor', 'start')
                .style('font-size', 10)
                .attr('dy', '1.2em')
                .attr('dx', '0.2em')
            .select(function() { return this.parentNode; })
                .append('clipPath')
                .attr('id', d => 'clip-'+d.data.id)
                    .append('use')
                    .attr('xlink:href', d => '#rect-'+d.data.id)
                .select(function() { return this.parentNode; })
            .select(function() { return this.parentNode; })
            .style('opacity', 0).transition(t_last).style('opacity', 1);

        catg.enter()
            .append('text')
            .attr('class', 'catg')
            .attr('id', d => d)
            .text(d => d)
            .attr('x', metadata.x1+10)
            .attr('y', d => category_y[d])
            .attr('dy', '0.8em')
            .call(d3.drag().on('drag', drag_catg).on('end', drag_catg_end))
            .style('opacity', 0).transition(t_last).style('opacity', 1);


        gene.enter()
            .filter(d => d.height === 0)  // leaves
            .append('text')
            .attr('class', 'gene')
            .attr('id', d => d.data.id)
            .attr('x', y_axis_leaves_x)
            .attr('y', d => y[d.data.id])
            .text(d => d.data.name)
            .attr('font-family', 'sans-serif')
            .style('font-weight', 300)
            .style('cursor', 'pointer')
            .style('text-anchor', (y_axis_leaves_position === 'before' ? 'end' : 'start'))
            .attr('dy', rect_height)
            .on('click', (d) => GeneCards(d.data.name))
            .call(d3.drag().on('drag', drag_gene).on('end', drag_gene_end))
            .style('opacity', 0).transition(t_last).style('opacity', 1);

        gene.enter()
            .filter(d => d.depth === 1)  // gene sets
            .append('g')
            .attr('class', 'gene')
            .attr('id', d => d.data.id)
            .attr('transform', d => 'translate('+y_axis_nodes_x+','+d.x1+')rotate(-90)')
            .call(d3.drag().on('drag', drag_gene).on('end', drag_gene_end))
                .append('rect')
                .attr('class', 'gene_set_box')
                .attr('id', d => 'rect-' + d.data.id)
                .attr('width', d => d.x1 - d.x0)
                .attr('height', d => d.y1 - d.y0)
                .style('fill-opacity', 0)
                .style('stroke', 'black')
                .style('stroke-dasharray', d => (y_axis_nodes_position === 'before' ? pointing_right(d) : pointing_left(d)))
            .select(function() { return this.parentNode; })
                .append('text')
                .attr('class', 'gene_set_label')
                .attr('clip-path', d => 'url(#clip-' + d.data.id + ')')
                .text(d => d.data.name)
                .attr('font-family', 'sans-serif')
                .style('text-anchor', 'start')
                .style('font-size', 10)
                .attr('dy', '1.2em')
                .attr('dx', '0.2em')
            .select(function() { return this.parentNode; })
                .append('clipPath')
                .attr('id', d => 'clip-'+d.data.id)
                    .append('use')
                    .attr('xlink:href', d => '#rect-'+d.data.id)
                .select(function() { return this.parentNode; })
            .select(function() { return this.parentNode; })
            .style('opacity', 0).transition(t_last).style('opacity', 1);



        rect.enter()
            .append('rect')
            .attr('class', 'rect')
            .attr('id', d => d.id)
            .attr('x', d => x[d.sample_id])
            .attr('y', d => y[d.gene_id])
            .attr('width', rect_width-2)
            .attr('height', rect_height-2)
            // .attr('fill', d => colors[values](d[values]))
            .style('opacity', 0).transition(t_last).style('opacity', 1);


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

        g.selectAll('.rect')
            .style('fill', (d) => colors[values](d[values]));



        if (show_legends) {

        } else {

        }

    }


    /////////////////////////////////////////////////////////////////////////////
                          ///////    Drag Axes    ///////
    /////////////////////////////////////////////////////////////////////////////

    function drag_node(d, hierarchy, xy, attr, svg_class) {

        let initial_position = (node) => (xy === 'x' ? node.x0 : (node.height === 0 ? node.x0 : node.x1));

        index_of_dragging_node = d.data.order;
        current_position_of_dragging_node = (d.x ? d.x : initial_position(d)) + (xy === 'x' ? d3.event.dx : d3.event.dy);
        dragging_node_width = d.x1 - d.x0;

        set_nodes = hierarchy.descendants().filter(node => node.parent && node.parent.data.id === d.parent.data.id);

        let expr = (node) => {

            index_of_other_node = node.data.order;
            original_position_of_other_node = initial_position(node);
            other_node_width = node.x1 - node.x0;

            if (index_of_other_node < index_of_dragging_node && current_position_of_dragging_node - (xy === 'y' ? dragging_node_width : other_node_width) < original_position_of_other_node) {
                return original_position_of_other_node + dragging_node_width + (margins[xy][d.height] || 0) + 2;
            }

            if (index_of_other_node > index_of_dragging_node && current_position_of_dragging_node > original_position_of_other_node - (xy === 'y' ? other_node_width : dragging_node_width)) {
                return original_position_of_other_node - dragging_node_width - (margins[xy][d.height] || 0) - 2;
            }

            if (index_of_other_node !== index_of_dragging_node) { return original_position_of_other_node; }


            if (index_of_other_node === index_of_dragging_node) {
                if (xy === 'y' && node.height > 0) { return clamp(node.parent.x0+other_node_width, node.parent.x1)(current_position_of_dragging_node); }
                else { return clamp(node.parent.x0, node.parent.x1-other_node_width)(current_position_of_dragging_node); }
            }
        };

        updated_xy = _.object(
            flatten(
                set_nodes.map(node => {
                    delta = expr(node) - initial_position(node);
                    return node.descendants().map(des => [des.data.id, initial_position(des)+delta]);
                })
            )
        );

        g.selectAll(svg_class).filter(node => node.height === 0 && node.data.id in updated_xy).each(node => {node.x = updated_xy[node.data.id]}).attr(xy, node => node.x);
        g.selectAll('.rect').filter(rect => updated_xy[rect[attr+'_id']]).attr(xy, (rect) => updated_xy[rect[attr+'_id']]);

        if (xy === 'x') {
            g.selectAll(svg_class).filter(node => node.height === 0 && node.data.id in updated_xy).attr('transform', node => 'rotate('+x_axis_leaves_rotation+','+node.x+','+x_axis_leaves_y+')');
        }

        if (d.height > 0) {
            g.selectAll(svg_class).filter(node => node.height > 0 && node.data.id in updated_xy).each(node => {node.x = updated_xy[node.data.id]}).attr('transform', function(node) {
                if (xy === 'x') { return 'translate('+node.x+','+category_y[node.data.category]+')'; }
                if (xy === 'y') { return 'translate('+y_axis_nodes_x+','+node.x+')rotate(-90)'; }
            });
        }

    }

    function drag_node_end(d, hierarchy, xy, ordered_this_ids, this_wise, other_wise) {

        set_nodes = hierarchy.descendants().filter(node => node.parent && node.parent.data.id === d.parent.data.id);

        new_order = _.object(set_nodes.filter(node => node.data.id !== d.data.id)
                                       .map(node => [node.x, node.data.id])
                                       .concat([[d.x, d.data.id]])
                                       .sort((a, b) => a[0] - b[0])
                                       .map(([y, id], i) => [id, i]));

        hierarchy.each(node => { if (node.data.id && node.data.id in new_order) { node.data.order = new_order[node.data.id]} });
        hierarchy.each(node => { node.x = undefined; }); // do I even need this?

        old_index = ordered_this_ids.indexOf(d.data.id);
        hierarchy.sort(function(a, b) { return b.height - a.height || a.data.order - b.data.order; });
        ordered_this_ids = hierarchy.leaves().map(node => node.data.id);   //// THIS IS A BUG -- IT ONLY CHANGES THE LOCAL COPY!!!
        new_index = ordered_this_ids.indexOf(d.data.id);

        this_wise.move(old_index, new_index);
        other_wise.forEach((other) => other.move(old_index, new_index));

        render();

    }

    // Drag category

    function drag_catg(d) {

        bounds = d3.extent(Object.values(category_y));
        current_position_of_dragging_category = clamp(bounds[0], bounds[1])(d3.event.y);
        depth = _.object(categories.map((catg, i) => [catg, i]));

        let expr = (category) => {

            if (depth[category] > depth[d] && current_position_of_dragging_category > category_y[category] - rect_height/2) { return category_y[category] - rect_height * (x_axis_nodes_position === 'before' ? 1 : -1); }

            if (depth[category] < depth[d] && current_position_of_dragging_category < category_y[category] + rect_height/2) { return category_y[category] + rect_height * (x_axis_nodes_position === 'before' ? 1 : -1); }

            if (category !== d) { return category_y[category]; }

            if (category === d) { return current_position_of_dragging_category; }

        }

        updated_category_y = _.object(categories, categories.map(expr));

        g.selectAll('.meta').filter(node => node.data.category in updated_category_y).each(node => {node.y = updated_category_y[node.data.category]}).attr('transform', node => 'translate('+node.x0+','+node.y+')');
        g.selectAll('.catg').attr('y', d => updated_category_y[d]);

    }

    function drag_catg_end(d) {

        category_y = {};
        g.selectAll('.catg').each(function(d) { category_y[d] = d3.select(this).attr('y'); });
        categories = Object.entries(category_y).sort((a, b) => a[1] - b[1]).map(([category, pos]) => category);
        order();
        render();

    }


    /////////////////////////////////////////////////////////////////////////////
                          ///////    Brush Axes    ///////
    /////////////////////////////////////////////////////////////////////////////


    function setFocus(d) {
    }

    function removeFocus(d) {
    }

    function GeneCards(d) { window.open('http://www.genecards.org/cgi-bin/carddisp.pl?gene='+d,'_blank') }


    /////////////////////////////////////////////////////////////////////////////
                          ///////   Zoom & Resize    ///////
    /////////////////////////////////////////////////////////////////////////////

    svg.call(d3.zoom().on('zoom', zoomed)).on('wheel.zoom', wheeled);

    transform = d3.zoomTransform(g);
    transform.x += margin.left;
    transform.y += margin.top;
    g.attr('transform', transform);

    function zoomed() {
        current_transform = d3.zoomTransform(g);
        current_transform.x += d3.event.sourceEvent.movementX;
        current_transform.y += d3.event.sourceEvent.movementY;
        g.attr('transform', current_transform);
    }

    function wheeled() {
        current_transform = d3.zoomTransform(g);
        if (d3.event.ctrlKey) {
            current_transform.k = clamp(0.1, 5)(current_transform.k - d3.event.deltaY * 0.01);
        } else {
            current_transform.y = clamp(-(ordered_gene_ids.length*rect_height-100), h)(current_transform.y - d3.event.deltaY);
        }
        g.attr('transform', current_transform);
    }

    function resize() {
        svg.attr('width', $('#graph-container').innerWidth()).attr('height', $('#graph-container').innerHeight());
        w = $('#graph-container').innerWidth() - (margin.left + margin.right);
        h = $('#graph-container').innerHeight() - (margin.top + margin.bottom);
    }

    d3.select(window).on('resize', resize)

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





