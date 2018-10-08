

let range = n => [...Array(n).keys()];  // underscore has a range()
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

    /////////////////////////////////////////////////////////////////////////////
                    ///////    Styling Variables    ///////
    /////////////////////////////////////////////////////////////////////////////

    var negative_color = '#0000cc';
    var middle_color = '#c0c0c0';
    var positive_color = '#cc0000';
    var show_legends = false;

    var category_colors = _.object(categories.map((category) => [category, d3.scaleOrdinal(d3.schemeCategory10)]))

    var colors = {};

    var spacing = 1;
    var rect_width = 16;
    var rect_height = 16;
    var max_font_size = 16;
    var margins = {'sample_id': {'2': 10}, 'gene_id': {'1': 10}};
    var styles = {
        'nodes': {
            'sample_id': {
                'fill': d => category_colors[d.data.category](d.data.name)},
            'gene_id': {
                'fill-opacity': 0,
                'stroke': 'black',
                'stroke-dasharray': d => (x_axis_nodes_position === 'before' ? pointing_right(d) : pointing_left(d))}
        },
        'leaves': { 'sample_id': {}, 'gene_id': {}, }
    };
    var text_styles = {
        'font-family': 'sans-serif',
        'font-weight': 300,
        'cursor': 'pointer',
        'text-anchor': 'start',
    };


    var y_axis_leaves_position = 'before';
    var y_axis_nodes_position = 'before';
    var y_axis_nodes_x_width = 16;
    var x_axis_leaves_position = 'after';
    var x_axis_nodes_position = 'before';
    var x_axis_nodes_y_height = 16;
    var y_axis_style = 'genes';
    var x_axis_style = 'metadata';
    var t = false;
    var x_categories = categories;
    var y_categories = ['Gene Set'];
    var y_axis_leaves_x, y_axis_nodes_x, x_axis_nodes_y, x_axis_leaves_y;
    var x_axis_leaves_rotation = (x_axis_leaves_position === 'before') ? -60 : 60;

    let text_max_width = (tree, font_size) => d3.max(tree.leaves().map(leaf => leaf.data.name.length)) * font_size;

    // position[display_style][nodes_or_leaves]: (params) => int;
    var axis_position = {
        'genes': {
            leaves: (nodes_position, leaves_position, this_tree, other_tree, layer_width, text_width) =>
                (leaves_position === 'before') ? -10 : other_tree.x1 + 10,
            nodes:  (nodes_position, leaves_position, this_tree, other_tree, layer_width, text_width) =>
                (nodes_position === 'before') ? -(layer_width*this_tree.height)-10 - ((leaves_position === 'before') ? text_width : 0) : other_tree.x1 - layer_width + ((leaves_position === 'before') ? 10 : text_width),
        },
        'metadata': {
            leaves: (nodes_position, leaves_position, this_tree, other_tree, layer_width, text_width) =>
                (leaves_position === 'before') ? ((nodes_position === 'before') ? -(layer_width*(this_tree.height-1))-20 : -10) : other_tree.x1 + ((nodes_position === 'before') ? 10 : (layer_width*(this_tree.height-1))+20),
            nodes:  (nodes_position, leaves_position, this_tree, other_tree, layer_width, text_width) =>
                (nodes_position === 'before') ? -(layer_width*this_tree.height)-10 : other_tree.x1 - layer_width + 10,
        }
    };

    /////////////////////////////////////////////////////////////////////////////
                          ///////    Set Up Chart    ///////
    /////////////////////////////////////////////////////////////////////////////

    var svg = d3.select('#graph-container').append('svg').attr('xmlns', 'http://www.w3.org/2000/svg').attr('xmlns:xlink', 'http://www.w3.org/1999/xlink');
    var g = svg.append('g');
    svg.style('cursor', 'move');

    var selected_gene_sets = {}
    var genes = [];
    var metadata = {};
    var gene_wise = [];
    var gene_wise_indexer = {};
    var ordered_gene_wise = [];
    var sample_wise = [];
    var sample_wise_indexer = {};
    var sample_to_sample_id = {};
    var x, y, x_category_y, y_category_x;
    var x_tree, y_tree, x_attr, y_attr;

    var drag_y        = (d) => drag_node(     d, y_tree, 'y', y_attr);
    var drag_y_end    = (d) => drag_node_end( d, y_tree, 'y', sample_wise, ordered_gene_wise);

    var drag_x        = (d) => drag_node(     d, x_tree, 'x', x_attr);
    var drag_x_end    = (d) => drag_node_end( d, x_tree, 'x', ordered_gene_wise, sample_wise);

    var drag_ycat     = (d) => drag_catg(     d, 'y', y_categories, y_category_x, y_axis_nodes_x_width);
    var drag_ycat_end = (d) => drag_catg_end( d, 'y');

    var drag_xcat     = (d) => drag_catg(     d, 'x', x_categories, x_category_y, x_axis_nodes_y_height);
    var drag_xcat_end = (d) => drag_catg_end( d, 'x');


    var legends = g.append('g').attr('class', 'legends');

    var rect_resizer = g.append('circle')
                        .attr('class', 'resizer')
                        .attr('id', 'rect_resizer')
                        .attr('r', 20)
                        .style('cursor', 'crosshair')
                        .style('opacity', 0)
                        .call(d3.drag().on('drag', drag_rect_resizer).on('end', render));

    var xtre_resizer = g.append('circle')
                        .attr('class', 'resizer')
                        .attr('id', 'xtre_resizer')
                        .attr('r', 20)
                        .style('cursor', 'crosshair')
                        .style('opacity', 0)
                        .call(d3.drag().on('drag', drag_xtre_resizer).on('end', render));

    var ytre_resizer = g.append('circle')
                        .attr('class', 'resizer')
                        .attr('id', 'ytre_resizer')
                        .attr('r', 20)
                        .style('cursor', 'crosshair')
                        .style('opacity', 0)
                        .call(d3.drag().on('drag', drag_ytre_resizer).on('end', render));

    /////////////////////////////////////////////////////////////////////////////
                          ///////    Methods    ///////
    /////////////////////////////////////////////////////////////////////////////

    function restart({selected_gene_sets_=selected_gene_sets}={}) {

        selected_gene_sets = selected_gene_sets_;
        selected_gene_sets = _.uniq(selected_gene_sets, false, (gs,i) => gs.gene_set_name ? gs.gene_set_name : i);  // gross, but how else do you keep all the nulls?
        keys = Object.keys(Object.values(samples_by_genes_matrix)[0]);  // genes included in matrix

        genes = d3.hierarchy({
                    'id': 'genes',
                    'children': selected_gene_sets.filter(gs => gs.gene_set_name !== null || keys.includes(gs.genes[0])).map((gs, i) => {
                        if (gs.gene_set_name === null) {
                            return {'gene_set':null, 'name':gs.genes[0], 'id':'other'+'_'+gs.genes[0]}
                        } else {
                            return {
                                'id': safeStr(gs.gene_set_name),
                                'name': gs.gene_set_name,
                                'category': 'Gene Set',
                                'children': _.uniq(gs.genes).filter(gene => keys.includes(gene)).map(gene => {
                                     return {'gene_set':gs.gene_set_name, 'name':gene, 'id':(safeStr(gs.gene_set_name) || 'other')+'_'+gene}
                                })
                            }}
                        })
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

    function order({
        values_=values,
        sorting_=sorting,
        minimum_nonzero_=minimum_nonzero}={}) {

        values = values_;
        sorting = sorting_;
        minimum_nonzero = minimum_nonzero_;
        value_accessor = value_accessors[values];

        // METADATA  // to order the levels: change the order of categories

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
        ordered_gene_wise = genes.leaves().map(leaf => gene_wise[gene_wise_indexer[leaf.data.name]].map(sample => Object.assign(sample, {'gene_id':leaf.data.id})));

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

                    if (set_gene_values.length > 1) {

                        if (sorting === 'complete') { permutation_order = reorder.optimal_leaf_order()(set_gene_values); } // get dendogram out?
                        else if (sorting === 'pc1') { permutation_order = reorder.sort_order(genes_pc1); }

                        set_gene_ids = reorder.stablepermute(set_gene_ids, permutation_order);
                        permutation_order = _.object(set_gene_ids.map((id, i) => [id, i]));
                        node.children.forEach((gene) => gene.data.order = permutation_order[gene.data.id]);

                    } else {
                        node.children[0].data.order = 0;
                    }
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

        metadata.count().sort(function(a, b) { return b.height - a.height || a.data.order - b.data.order; });

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

    function set_transposition(t_) {

        t = t_ ? true : false;

        // [y_axis_style, x_axis_style] = [x_axis_style, y_axis_style];

        if (t) {

            x_categories = (genes.height > 1 ? ['Gene Set'] : []);
            y_categories = categories;

            x_tree = genes;
            x_attr = 'gene_id';
            y_tree = metadata;
            y_attr = 'sample_id';

            drag_y_end = (d) => drag_node_end( d, y_tree, 'y', sample_wise, ordered_gene_wise);
            drag_x_end = (d) => drag_node_end( d, x_tree, 'x', ordered_gene_wise, sample_wise);

        } else {

            x_categories = categories;
            y_categories = (genes.height > 1 ? ['Gene Set'] : []);

            x_tree = metadata;
            x_attr = 'sample_id';
            y_tree = genes;
            y_attr = 'gene_id';

            drag_y_end = (d) => drag_node_end( d, y_tree, 'y', ordered_gene_wise, sample_wise);
            drag_x_end = (d) => drag_node_end( d, x_tree, 'x', sample_wise, ordered_gene_wise);

        }
    }

    function position() {

        set_transposition(t);

        x_tree_across = (rect_width*x_tree.leaves().length)+spacing;
        x_tree_topdown = x_axis_nodes_y_height*(x_tree.height+1);
        d3.partition().size([x_tree_across, x_tree_topdown]).padding(spacing)(x_tree);
        offset(x_tree, margins[x_attr]);

        y_tree_across = (rect_height*y_tree.leaves().length)+spacing;
        y_tree_topdown = y_axis_nodes_x_width*(y_tree.height+1);
        d3.partition().size([y_tree_across, y_tree_topdown]).padding(spacing)(y_tree);
        offset(y_tree, margins[y_attr]);

        x = _.object(x_tree.leaves().map(leaf => [leaf.data.id, leaf.x0]));
        y = _.object(y_tree.leaves().map(leaf => [leaf.data.id, leaf.x0]));

        y_font_size = Math.min(rect_height-spacing, max_font_size);
        x_font_size = Math.min(rect_width-spacing, max_font_size);

        y_axis_leaves_x = axis_position[y_axis_style]['leaves'](y_axis_nodes_position, y_axis_leaves_position, y_tree, x_tree, y_axis_nodes_x_width, text_max_width(y_tree, y_font_size));
        y_axis_nodes_x  = axis_position[y_axis_style]['nodes']( y_axis_nodes_position, y_axis_leaves_position, y_tree, x_tree, y_axis_nodes_x_width, text_max_width(y_tree, y_font_size));
        x_axis_leaves_y = axis_position[x_axis_style]['leaves'](x_axis_nodes_position, x_axis_leaves_position, x_tree, y_tree, x_axis_nodes_y_height, text_max_width(x_tree, x_font_size)*2/3);
        x_axis_nodes_y  = axis_position[x_axis_style]['nodes']( x_axis_nodes_position, x_axis_leaves_position, x_tree, y_tree, x_axis_nodes_y_height, text_max_width(x_tree, x_font_size)*2/3);

        x_axis_leaves_rotation = (x_axis_leaves_position === 'before') ? -60 : 60;

        x_category_y = _.object(x_categories.map((c, i) => [c, x_axis_nodes_y + x_axis_nodes_y_height*(i+1)]));
        y_category_x = _.object(y_categories.map((c, i) => [c, y_axis_nodes_x + y_axis_nodes_x_width*(i+1)]));

    }

    function render({
        y_axis_leaves_position_=y_axis_leaves_position,
        y_axis_nodes_position_=y_axis_nodes_position,
        y_axis_style_=y_axis_style,
        x_axis_leaves_position_=x_axis_leaves_position,
        x_axis_nodes_position_=x_axis_nodes_position,
        x_axis_style_=x_axis_style,}={}) {

        y_axis_leaves_position = y_axis_leaves_position_;
        y_axis_nodes_position = y_axis_nodes_position_;
        y_axis_style = y_axis_style_;
        x_axis_leaves_position = x_axis_leaves_position_;
        x_axis_nodes_position = x_axis_nodes_position_;
        x_axis_style = x_axis_style_;

        position();

        rect = g.selectAll('.rect').data(flatten(ordered_gene_wise), d => d.id);
        ytre = g.selectAll('.ytre').data(y_tree.descendants(), d => d.data.id);
        xtre = g.selectAll('.xtre').data(x_tree.descendants(), d => d.data.id);
        xcat = g.selectAll('.xcat').data(x_categories, d => d);
        ycat = g.selectAll('.ycat').data(y_categories, d => d);

        // phase 1
            // rectangles which are exiting fade out
            // gene names which are exiting fade out
            // gene groups which are exiting fade out
        t_last = d3.transition().duration(500);
        if (rect.exit().size() > 0 || ytre.exit().size() > 0 || xtre.exit().size() > 0 || xcat.exit().size() > 0 || ycat.exit().size() > 0) {
            rect.exit().transition(t_last).style('opacity', 0).remove();
            ytre.exit().transition(t_last).style('opacity', 0).remove();
            xtre.exit().transition(t_last).style('opacity', 0).remove();
            ycat.exit().transition(t_last).style('opacity', 0).remove();
            xcat.exit().transition(t_last).style('opacity', 0).remove();
            t_last = t_last.transition().duration(500);
        }

        // phase 2
            // re-arrange ROWS
        rect.transition(t_last).attr('y', d => y[d[y_attr]]).attr('height', rect_height-spacing);
        ytre.filter(node => node.height === 0).transition(t_last).attr('y', d => y[d.data.id])
                                                                 .attr('x', y_axis_leaves_x)
                                                                 .style('text-anchor', (y_axis_leaves_position === 'before' ? 'end' : 'start'))
                                                                 .style('font-size', y_font_size)
                                                                 .attr('dy', y_font_size);
        ytre.filter(node => node.depth > 0 && node.height > 0).transition(t_last).attr('transform', d => 'translate('+(y_axis_nodes_x + d.y0)+','+d.x1+')rotate(-90)');
        ytre.filter(node => node.depth > 0 && node.height > 0).select('.ytre_box').transition(t_last).attr('width', d => d.x1 - d.x0).attr('height', d => d.y1 - d.y0).style('stroke-dasharray', d => (y_axis_nodes_position === 'before' ? pointing_right(d) : pointing_left(d)));
        ycat.transition(t_last).attr('y', x_axis_leaves_y).attr('x', d => y_category_x[d]).attr('dy', (x_axis_leaves_position === 'before' ? y_font_size : 0))
                               .attr('transform', d => 'rotate('+x_axis_leaves_rotation+','+y_category_x[d]+','+x_axis_leaves_y+')');
        t_last = t_last.transition().duration(500);

        // phase 3
            // re-arrange COLUMNS
        rect.transition(t_last).attr('x', d => x[d[x_attr]]).attr('width', rect_width-spacing);
        xtre.filter(node => node.height === 0).transition(t_last).attr('x', d => d.x0)
                                                                 .attr('y', x_axis_leaves_y)
                                                                 .attr('transform', d => 'rotate('+x_axis_leaves_rotation+','+d.x0+','+x_axis_leaves_y+')')
                                                                 .style('font-size', x_font_size)
                                                                 .attr('dy', (x_axis_leaves_position === 'before' ? x_font_size : 0));
        xtre.filter(node => node.depth > 0 && node.height > 0).transition(t_last).attr('transform', d => 'translate('+d.x0+','+(x_axis_nodes_y + d.y0)+')');
        xtre.filter(node => node.depth > 0 && node.height > 0).select('.xtre_box').transition(t_last).attr('width', d => d.x1 - d.x0).attr('height', d => d.y1 - d.y0).style('stroke-dasharray', d => (y_axis_nodes_position === 'before' ? pointing_right(d) : pointing_left(d)));
        xcat.transition(t_last).attr('x', y_axis_leaves_x).attr('y', d => x_category_y[d]).style('text-anchor', (y_axis_leaves_position === 'before' ? 'end' : 'start'));
        t_last = t_last.transition().duration(500);

        // phase 4
            // rectangles which are entering get appended
            // gene names / groups which are entering get appended
            // sample names / groups which are entering get appended
        xtre.enter()
            .filter(d => d.height === 0)  // leaves
            .append('text')
            .attr('class', 'xtre')
            .attr('id', d => d.data.id)
            .attr('x', d => x[d.data.id])
            .attr('y', x_axis_leaves_y)
            .attr('dy', (x_axis_leaves_position === 'before' ? x_font_size : 0))
            .attr('transform', d => 'rotate('+x_axis_leaves_rotation+','+d.x0+','+x_axis_leaves_y+')')
            .text(d => d.data.name)
            .styles(text_styles)
            .style('font-size', x_font_size)
            .call(d3.drag().on('drag', drag_x).on('end', drag_x_end))
            .style('opacity', 0).transition(t_last).style('opacity', 1);

        xtre.enter()
            .filter(d => d.children !== undefined && d.depth > 0)  // internal nodes
            .append('g')
            .attr('class', 'xtre')
            .attr('id', d => d.data.id)
            .attr('transform', d => 'translate('+d.x0+','+(x_axis_nodes_y + d.y0)+')')
            .call(d3.drag().on('drag', drag_x).on('end', drag_x_end))
                .append('rect')
                .attr('class', 'xtre_box')
                .attr('id', d => 'rect-' + d.data.id)
                .attr('width', d => d.x1 - d.x0)
                .attr('height', d => d.y1 - d.y0)
                .styles(styles['nodes'][x_attr])
            .select(function() { return this.parentNode; })
                .append('text')
                .attr('class', 'xtre_label')
                .attr('clip-path', d => 'url(#clip-' + d.data.id + ')')
                .text(d => d.data.name)
                .styles(text_styles)
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

        xcat.enter()
            .append('text')
            .attr('class', 'xcat')
            .attr('id', d => d)
            .text(d => d)
            .attr('x', y_axis_leaves_x)
            .attr('y', d => x_category_y[d])
            .attr('dy', '0.8em')
            .style('font-size', Math.min(x_axis_nodes_y_height-spacing, max_font_size))
            .styles(text_styles)
            .style('text-anchor', (y_axis_leaves_position === 'before' ? 'end' : 'start'))
            .call(d3.drag().on('drag', drag_xcat).on('end', drag_xcat_end))
            .style('opacity', 0).transition(t_last).style('opacity', 1);


        ytre.enter()
            .filter(d => d.height === 0)  // leaves
            .append('text')
            .attr('class', 'ytre')
            .attr('id', d => d.data.id)
            .attr('x', y_axis_leaves_x)
            .attr('y', d => y[d.data.id])
            .text(d => d.data.name)
            .style('font-size', y_font_size)
            .styles(text_styles)
            .style('text-anchor', (y_axis_leaves_position === 'before' ? 'end' : 'start'))
            .attr('dy', y_font_size)
            .on('click', (d) => GeneCards(d.data.name))
            .call(d3.drag().on('drag', drag_y).on('end', drag_y_end))
            .style('opacity', 0).transition(t_last).style('opacity', 1);

        ytre.enter()
            .filter(d => d.children !== undefined && d.depth > 0)  // gene sets
            .append('g')
            .attr('class', 'ytre')
            .attr('id', d => d.data.id)
            .attr('transform', d => 'translate('+(y_axis_nodes_x + d.y0)+','+d.x1+')rotate(-90)')
            .call(d3.drag().on('drag', drag_y).on('end', drag_y_end))
                .append('rect')
                .attr('class', 'ytre_box')
                .attr('id', d => 'rect-' + d.data.id)
                .attr('width', d => d.x1 - d.x0)
                .attr('height', d => d.y1 - d.y0)
                .styles(styles['nodes'][y_attr])
            .select(function() { return this.parentNode; })
                .append('text')
                .attr('class', 'ytre_label')
                .attr('clip-path', d => 'url(#clip-' + d.data.id + ')')
                .text(d => d.data.name)
                .styles(text_styles)
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

        ycat.enter()
            .append('text')
            .attr('class', 'ycat')
            .attr('id', d => d)
            .text(d => d)
            .attr('x', d => y_category_x[d])
            .attr('y', x_axis_leaves_y)
            .attr('transform', d => 'rotate('+x_axis_leaves_rotation+','+y_category_x[d]+','+x_axis_leaves_y+')')
            .style('font-size', Math.min(x_axis_nodes_y_height-spacing, max_font_size))
            .styles(text_styles)
            .attr('dx', '0.2em')
            .call(d3.drag().on('drag', drag_ycat).on('end', drag_ycat_end))
            .style('opacity', 0).transition(t_last).style('opacity', 1);


        rect.enter()
            .append('rect')
            .attr('class', 'rect')
            .attr('id', d => d.id)
            .attr('x', d => x[d[x_attr]])
            .attr('y', d => y[d[y_attr]])
            .attr('width', rect_width-spacing)
            .attr('height', rect_height-spacing)
            // .attr('fill', d => colors[values](d[values]))
            .style('opacity', 0).transition(t_last).style('opacity', 1);


        g.select('#rect_resizer').attr('cx', d3.max(Object.values(x))+rect_width).attr('cy', d3.max(Object.values(y))+rect_height);
        g.select('#xtre_resizer').attr('cx', d3.max(Object.values(x))+rect_width).attr('cy', x_axis_nodes_y + (x_axis_nodes_position === 'before' ? x_axis_nodes_y_height : x_axis_nodes_y_height*x_tree.height));
        g.select('#ytre_resizer').attr('cx', y_axis_nodes_x + (y_axis_nodes_position === 'before' ? y_axis_nodes_x_width : y_axis_nodes_x_width*y_tree.height)).attr('cy', d3.max(Object.values(y))+rect_height);

        style();

    }

    function style({negative_color_=negative_color,
                    middle_color_=middle_color,
                    positive_color_=positive_color,
                    show_legends_=show_legends}={}) {

        negative_color = negative_color_,
        middle_color = middle_color_,
        positive_color = positive_color_,
        show_legends = show_legends_;

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

    function resize_fig({
        rect_width_=rect_width,
        rect_height_=rect_height,
        x_axis_nodes_y_height_=x_axis_nodes_y_height,
        y_axis_nodes_x_width_=y_axis_nodes_x_width,}={}) {

        rect_width = rect_width_;
        rect_height = rect_height_;
        x_axis_nodes_y_height = x_axis_nodes_y_height_;
        y_axis_nodes_x_width = y_axis_nodes_x_width_;

        position();

        rect = g.selectAll('.rect');
        ytre = g.selectAll('.ytre');
        xtre = g.selectAll('.xtre');
        ycat = g.selectAll('.ycat');
        xcat = g.selectAll('.xcat');

        rect.attr('y', d => y[d[y_attr]])
            .attr('x', d => x[d[x_attr]])
            .attr('width', rect_width-spacing)
            .attr('height', rect_height-spacing);

        ytre.filter(node => node.height === 0).attr('y', d => y[d.data.id])
                                              .attr('x', y_axis_leaves_x)
                                              .style('font-size', y_font_size)
                                              .attr('dy', y_font_size);
        ytre.filter(node => node.depth > 0 && node.height > 0).attr('transform', d => 'translate('+(y_axis_nodes_x + d.y0)+','+d.x1+')rotate(-90)')
        ytre.filter(node => node.depth > 0 && node.height > 0).select('.ytre_box').attr('width', d => d.x1 - d.x0).attr('height', d => d.y1 - d.y0).style('stroke-dasharray', d => (y_axis_nodes_position === 'before' ? pointing_right(d) : pointing_left(d)));
        ycat.attr('y', x_axis_leaves_y)
            .attr('x', d => y_category_x[d])
            .attr('transform', d => 'rotate('+x_axis_leaves_rotation+','+y_category_x[d]+','+x_axis_leaves_y+')');

        xtre.filter(node => node.height === 0).attr('x', d => d.x0)
                                              .attr('y', x_axis_leaves_y)
                                              .attr('transform', d => 'rotate('+x_axis_leaves_rotation+','+d.x0+','+x_axis_leaves_y+')')
                                              .style('font-size', x_font_size)
                                              .attr('dy', (x_axis_leaves_position === 'before' ? x_font_size : 0))
        xtre.filter(node => node.depth > 0 && node.height > 0).attr('transform', d => 'translate('+d.x0+','+(x_axis_nodes_y + d.y0)+')');
        xtre.filter(node => node.depth > 0 && node.height > 0).select('.xtre_box').attr('width', d => d.x1 - d.x0).attr('height', d => d.y1 - d.y0).style('stroke-dasharray', d => (y_axis_nodes_position === 'before' ? pointing_right(d) : pointing_left(d)));
        xcat.attr('x', y_axis_leaves_x)
            .attr('y', d => x_category_y[d]);

    }


    ///////////////////////////////////////////////////////////////////////////
                          ///////      Drag      ///////
    ///////////////////////////////////////////////////////////////////////////

    function drag_node(d, hierarchy, xy, attr) {

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
                return original_position_of_other_node + dragging_node_width + (margins[attr][d.height] || 0) + spacing;
            }

            if (index_of_other_node > index_of_dragging_node && current_position_of_dragging_node > original_position_of_other_node - (xy === 'y' ? other_node_width : dragging_node_width)) {
                return original_position_of_other_node - dragging_node_width - (margins[attr][d.height] || 0) - spacing;
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

        g.selectAll('.'+xy+'tre').filter(node => node.height === 0 && node.data.id in updated_xy).each(node => {node.x = updated_xy[node.data.id]}).attr(xy, node => node.x);
        g.selectAll('.rect').filter(rect => updated_xy[rect[attr]]).attr(xy, (rect) => updated_xy[rect[attr]]);

        if (xy === 'x') {
            g.selectAll('.'+xy+'tre').filter(node => node.height === 0 && node.data.id in updated_xy).attr('transform', node => 'rotate('+x_axis_leaves_rotation+','+node.x+','+x_axis_leaves_y+')');
        }

        if (d.height > 0) {
            g.selectAll('.'+xy+'tre').filter(node => node.height > 0 && node.data.id in updated_xy).each(node => {node.x = updated_xy[node.data.id]}).attr('transform', function(node) {
                if (xy === 'x') { return 'translate('+node.x+','+(x_axis_nodes_y + node.y0)+')'; }
                if (xy === 'y') { return 'translate('+(y_axis_nodes_x+node.y0)+','+node.x+')rotate(-90)'; }
            });
        }

    }

    function drag_node_end(d, hierarchy, xy, this_wise, other_wise) {

        set_nodes = hierarchy.descendants().filter(node => node.parent && node.parent.data.id === d.parent.data.id);

        new_order = _.object(set_nodes.filter(node => node.data.id !== d.data.id)
                                       .map(node => [node.x, node.data.id])
                                       .concat([[d.x, d.data.id]])
                                       .sort((a, b) => a[0] - b[0])
                                       .map(([y, id], i) => [id, i]));

        hierarchy.each(node => { if (node.data.id && node.data.id in new_order) { node.data.order = new_order[node.data.id]} });
        hierarchy.each(node => { node.x = undefined; }); // do I even need this?

        old_index = hierarchy.leaves().map(leaf => leaf.data.id).indexOf(d.data.id);
        hierarchy.sort(function(a, b) { return b.height - a.height || a.data.order - b.data.order; });
        new_index = hierarchy.leaves().map(leaf => leaf.data.id).indexOf(d.data.id);

        this_wise.move(old_index, new_index);
        other_wise.forEach((other) => other.move(old_index, new_index));

        render();

    }

    // Drag category

    function drag_catg(d, xy, xy_categories, xy_categories_y, layer_width) {

        yx = (xy === 'x' ? 'y' : 'x');
        bounds = d3.extent(Object.values(xy_categories_y));
        current_position_of_dragging_category = clamp(bounds[0], bounds[1])(d3.event[yx]);
        depth = _.object(xy_categories.map((catg, i) => [catg, i]));

        let expr = (category) => {

            if (depth[category] > depth[d] && current_position_of_dragging_category > xy_categories_y[category] - layer_width/2) { return xy_categories_y[category] - layer_width; }

            if (depth[category] < depth[d] && current_position_of_dragging_category < xy_categories_y[category] + layer_width/2) { return xy_categories_y[category] + layer_width; }

            if (category !== d) { return xy_categories_y[category]; }

            if (category === d) { return current_position_of_dragging_category; }

        }

        updated_xy_categories_y = _.object(xy_categories, xy_categories.map(expr));

        g.selectAll('.'+xy+'cat').attr(yx, d => updated_xy_categories_y[d]);

        if (xy === 'y') {
            g.selectAll('.ycat').attr('transform', d => 'rotate('+x_axis_leaves_rotation+','+updated_xy_categories_y[d]+','+x_axis_leaves_y+')')
        }

        g.selectAll('.'+xy+'tre').filter(node => node.data.category in updated_xy_categories_y).each(node => {node.y = updated_xy_categories_y[node.data.category]}).attr('transform', function(node) {
            if (xy === 'x') { return 'translate('+node.x0+','+node.y+')'; }
            if (xy === 'y') { return 'translate('+node.y+','+node.x1+')rotate(-90)'; }
        });

    }

    function drag_catg_end(d, xy) {

        yx = (xy === 'x' ? 'y' : 'x');
        xy_categories_y = {};
        g.selectAll('.'+xy+'cat').each(function(d) { xy_categories_y[d] = d3.select(this).attr(yx); });
        updated_categories = Object.entries(xy_categories_y).sort((a, b) => a[1] - b[1]).map(([category, pos]) => category);
        if (_.isEqual(_.sortBy(categories), _.sortBy(updated_categories))) { categories = updated_categories; }
        order();
        render();

    }

    // Drag Resizers

    function drag_rect_resizer(d) {

        d3.select(this).attr('cx', d3.event.x).attr('cy', d3.event.y);
        resize_fig({
            'rect_width_': Math.max((d3.event.x - x_tree.leaves().last().offset), x_tree.leaves().length*2) / x_tree.leaves().length,
            'rect_height_': Math.max((d3.event.y - y_tree.leaves().last().offset), y_tree.leaves().length*2) / y_tree.leaves().length,
        });
    }

    function drag_xtre_resizer(d) {

        d3.select(this).attr('cx', d3.event.x).attr('cy', d3.event.y);
        resize_fig({
            'x_axis_nodes_y_height_': x_axis_nodes_y_height + (x_axis_nodes_position === 'before' ? -1 : 1) * (d3.event.dy / (x_tree.height-1))
        });
    }

    function drag_ytre_resizer(d) {

        d3.select(this).attr('cx', d3.event.x).attr('cy', d3.event.y);
        resize_fig({
            'y_axis_nodes_x_width_': y_axis_nodes_x_width + (y_axis_nodes_position === 'before' ? -1 : 1) * (d3.event.dx / (y_tree.height-1))
        });
    }



    /////////////////////////////////////////////////////////////////////////////
                          ///////    Hover    ///////
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
            current_transform.y = clamp(-(genes.leaves().length*rect_height-100), h)(current_transform.y - d3.event.deltaY);
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

        transpose     : function() { t = !t; [rect_width, rect_height] = [rect_height, rect_width]; render(); },

        get_sorted_genes: () => genes,

        set_reordering: (reordering_) => { reordering = reordering_; if (reordering) { order(); render(); } },
    }

}





