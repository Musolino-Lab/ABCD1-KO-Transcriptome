

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
Array.prototype.insert = function(index, item) { this.splice( index, 0, item ); return this; };

function Heatmap(samples_by_genes_matrix, gene_sets, classes, separate_zscore_by, refresh_genes_cb) {

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

    var samples_to_bin = _(classes).mapObject(categories_to_values => Object.entries(_(categories_to_values).pick(separate_zscore_by)).sort().reduce((acc, [category, value]) => (acc ? acc+'-'+value : value), ''));
    var bin_to_samples = _(Object.keys(samples_to_bin)).groupBy(sample => samples_to_bin[sample]);

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


    /////////////////////////////////////////////////////////////////////////////
                    ///////    Styling Variables    ///////
    /////////////////////////////////////////////////////////////////////////////

    var negative_color = '#0000cc';
    var middle_color = '#c0c0c0';
    var positive_color = '#cc0000';
    var show_legends = false;

    var category_colors = _.object(categories.map((category) => [category, d3.scaleOrdinal(d3.schemeCategory10)]))

    var color_style = 'interpolateTriplet'
    var colors = null;

    var spacing = 1;
    var rect_width = 16;
    var rect_height = 16;

    var margins = {'sample_id': {'2': 10}, 'gene_id': {'1': 10}};
    var styles = {
        'nodes': {
            'sample_id': {
                'fill': d => category_colors[d.data.category](d.data.name)},
            'gene_id': {
                'fill-opacity': 0,
                'stroke': 'black',
                'stroke-dasharray': d => pointing_right(d)}
        },
        'leaves': { 'sample_id': {}, 'gene_id': {}, }
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
    var rotation = 60;
    var x_axis_leaves_rotation = (x_axis_leaves_position === 'before') ? -rotation : rotation;
    var y_font_size, x_font_size, x_cat_font_size, y_cat_font_size, xtre_label_font_size, ytre_label_font_size;
    var show_x_level_names = true, show_y_level_names = true;

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

    var max_font_size = 16;
    var text_styles = {
        'font-family': 'sans-serif',
        'font-weight': 300,
        'cursor': 'pointer',
        'text-anchor': 'start',
    };

    let text_max_width = (tree, font_size) => (d3.max(tree.leaves().filter(leaf => leaf.depth > 0).map(leaf => leaf.data.name.length)) || 0) * font_size;

    var clear_styles = {
        'fill': 'red',
        'stroke': 'black',
        'cursor': 'pointer',
        'opacity': 0,
    }

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
    // var sample_wise_indexer = {};
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
        console.log('restart', selected_gene_sets);
        selected_gene_sets = _.uniq(selected_gene_sets, false, (gs,i) => gs.gene_set_name ? gs.gene_set_name : i);  // gross, but how else do you keep all the nulls?
        keys = Object.keys(Object.values(samples_by_genes_matrix)[0]);  // genes included in matrix

        genes = d3.hierarchy({
                    'id': 'genes',
                    'children': selected_gene_sets.filter(gs => _.any(gs.genes, gene => keys.includes(gene))).map((gs, i) => {
                        if (gs.gene_set_name === null) {
                            return {'gene_set':null, 'name':gs.genes[0], 'order':i, 'id':'other'+'_'+gs.genes[0]}
                        } else {
                            return {
                                'id': safeStr(gs.gene_set_name),
                                'name': gs.gene_set_name,
                                'order':i,
                                'category': 'Gene Set',
                                'children': _.uniq(gs.genes).filter(gene => keys.includes(gene)).map((gene, i) => {
                                     return {'gene_set':gs.gene_set_name, 'name':gene, 'order':i, 'id':safeStr(gs.gene_set_name)+'_'+gene}
                                })
                            }}
                        })
                });
        console.log(genes);

        if (genes.data.children.length === 0) { return clear_fig(); }

        matrix = _(samples_by_genes_matrix).mapObject((sample) => _(sample).pick(genes.leaves().map(d => d.data.name)));

        sample_wise = Object.entries(matrix).map(([sample, genes]) =>
            Object.entries(genes).map(([gene, count]) => { return {
                'id'        : sample+'_'+gene,
                'sample'    : sample,
                'sample_id' : null,  // assigned later
                'gene'      : gene,
                'gene_id'   : null,  // assigned later
                'count'     : count,
                'logcount'  : Math.log10(count+1),
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
        // sample_wise_indexer = _.object(sample_wise.map((sample, i) => [sample[0].sample, i]));

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

        // METADATA  // to order the levels: change the order of categories

        hierarchy = {'id':'metadata', 'children':[]};
        Object.entries(classes).filter(([sample_id, metadata]) => sample_id in samples_by_genes_matrix).forEach(([sample_id, metadata]) => {
            pointer = hierarchy.children;
            prefix = 'metadata';
            categories.forEach((category, depth) => {
                value = metadata[category];
                prefix += '-'+safeStr(value);
                existing_index_for_value = _(pointer).findIndex({'id':prefix});
                if (existing_index_for_value > -1) {
                    if (depth+1 === categories.length) { pointer[existing_index_for_value].children.push({'id':prefix+'-'+sample_id, 'name':sample_id, 'order':pointer[existing_index_for_value].children.length}); }
                    else { pointer = pointer[existing_index_for_value].children; }
                } else {
                    if (depth+1 === categories.length) { pointer.push({'id':prefix,'name':value,'category':category,'order':pointer.length,'children':[{'id':prefix+'-'+sample_id, 'name':sample_id, 'order':0}]}); }
                    else {
                        pointer.push({'id':prefix,'name':value,'category':category,'order':pointer.length,'children':[]});
                        pointer = pointer[_(pointer).findIndex({'id':prefix})].children;
                    }
                }
            })
        });

        metadata = d3.hierarchy(hierarchy);

        sample_to_sample_id = _.object(metadata.leaves().map(leaf => [leaf.data.name, leaf.data.id]));
        sample_wise.forEach(by_sample => by_sample.forEach(sample => sample.sample_id = sample_to_sample_id[sample.sample]));

        // GENES

        if (gene_wise.length === 0) { return; }

        else {
            // Filter by number non-zeros
            ordered_gene_wise = genes.leaves().map(leaf => gene_wise[gene_wise_indexer[leaf.data.name]].map(sample => Object.assign(sample, {'gene_id':leaf.data.id})));

            ordered_gene_wise = ordered_gene_wise.filter((gene) => gene[0].num_nonzeros >= minimum_nonzero);

            if (ordered_gene_wise.length === 0) { return; }  // do something smart here.
        }

        if (reordering && ordered_gene_wise.length > 1) {

            reorder_leaves(genes, ordered_gene_wise, 'gene_id');
            reorder_leaves(metadata, sample_wise, 'sample_id');

        }

        genes.each(node => console.log(node.data.id, node.data.order));

        genes.count().sort(function(a, b) { return a.depth - b.depth || a.data.order - b.data.order; });

        metadata.count().sort(function(a, b) { return a.depth - b.depth || a.data.order - b.data.order; });

    }

    function reorder_leaves(hierarchy, this_wise, id_attr) {

        hierarchy.each(node => {
            if (node.height === 1) {  // for each gene set / smallest set of samples

                set_leaf_ids = node.children.map(leaf => leaf.data.id);
                set_leaves = this_wise.filter(vec => set_leaf_ids.includes(vec[0][id_attr]));

                if (set_leaves.length > 1) {

                    set_leaf_values = set_leaves.map(vec => vec.map(leaf => leaf[values] || 0));

                    // Set leaf-Wise PC1
                    leaves_pc1 = reorder.pca1d(reorder.transpose(set_leaf_values));
                    _.zip(set_leaves, leaves_pc1).forEach(([vec, pc1]) => vec.forEach(d => d.pc1 = pc1));

                    if (sorting === 'complete') { permutation_order = reorder.optimal_leaf_order()(set_leaf_values); } // get dendogram out?
                    else if (sorting === 'pc1') { permutation_order = reorder.sort_order(leaves_pc1); }

                    set_leaf_ids = reorder.stablepermute(set_leaf_ids, permutation_order);
                    permutation_order = _.object(set_leaf_ids.map((id, i) => [id, i]));
                    node.children.forEach((leaf) => leaf.data.order = permutation_order[leaf.data.id]);

                } else {
                    node.children[0].data.order = 0;
                }
            }
        });

    }

    function offset(partition, margin) {

        counter = 0; current_depth = 0;
        partition.each((node) => {
            node.num_below = _(node.descendants().slice(1).map(d => d.height)).countBy();

            if (node.depth !== current_depth) { current_depth = node.depth; counter = 0; }
            node.sibling_index = counter; counter += 1;
        });

        partition.offset = 0;
        partition.num_left = 0;

        partition.each((node) => {
            if (node.parent) {
                left_sibling = _(node.parent.children.filter(sibling => sibling.sibling_index < node.sibling_index)).sortBy(sibling => sibling.sibling_index).last();

                if (!left_sibling) {
                    node.offset = node.parent.offset;
                } else {
                    node.offset = (
                        left_sibling.offset +
                        d3.sum(Object.entries(left_sibling.num_below).map(([level, num]) => (margin[level] || 0)*(num-1) )) +
                        d3.max([(margin[node.height] || 0), (margin[left_sibling.height] || 0)])
                    );
                }

                node.x0 += node.offset;
                node.x1 += node.offset + d3.sum(Object.entries(node.num_below).map(([level, num]) => (margin[level] || 0)*(num-1)));
            }
        });

        partition.x1 += d3.max(partition.descendants().map(node => node.offset));

    }

    function set_transposition(t_) {

        t = t_ ? true : false;

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
        x_cat_font_size = Math.min(x_axis_nodes_y_height-spacing, max_font_size);
        y_cat_font_size = Math.min(y_axis_nodes_x_width-spacing, max_font_size);
        xtre_label_font_size = Math.min(x_axis_nodes_y_height*2/3-spacing, max_font_size);
        ytre_label_font_size = Math.min(y_axis_nodes_x_width*2/3-spacing, max_font_size);

        y_axis_leaves_x = axis_position[y_axis_style]['leaves'](y_axis_nodes_position, y_axis_leaves_position, y_tree, x_tree, y_axis_nodes_x_width, text_max_width(y_tree, y_font_size));
        y_axis_nodes_x  = axis_position[y_axis_style]['nodes']( y_axis_nodes_position, y_axis_leaves_position, y_tree, x_tree, y_axis_nodes_x_width, text_max_width(y_tree, y_font_size));
        x_axis_leaves_y = axis_position[x_axis_style]['leaves'](x_axis_nodes_position, x_axis_leaves_position, x_tree, y_tree, x_axis_nodes_y_height, text_max_width(x_tree, x_font_size)*2/3);
        x_axis_nodes_y  = axis_position[x_axis_style]['nodes']( x_axis_nodes_position, x_axis_leaves_position, x_tree, y_tree, x_axis_nodes_y_height, text_max_width(x_tree, x_font_size)*2/3);

        x_axis_leaves_rotation = (x_axis_leaves_position === 'before') ? -rotation : rotation;

        x_category_y = _.object(x_categories.map((c, i) => [c, x_axis_nodes_y + x_axis_nodes_y_height*(i+1)]));
        y_category_x = _.object(y_categories.map((c, i) => [c, y_axis_nodes_x + y_axis_nodes_x_width*(i+1)]));

    }

    function render({
        spacing_=spacing,
        y_axis_leaves_position_=y_axis_leaves_position,
        y_axis_nodes_position_=y_axis_nodes_position,
        y_axis_style_=y_axis_style,
        x_axis_leaves_position_=x_axis_leaves_position,
        x_axis_nodes_position_=x_axis_nodes_position,
        x_axis_style_=x_axis_style,}={}) {

        spacing = spacing_;
        y_axis_leaves_position = y_axis_leaves_position_;
        y_axis_nodes_position = y_axis_nodes_position_;
        y_axis_style = y_axis_style_;
        x_axis_leaves_position = x_axis_leaves_position_;
        x_axis_nodes_position = x_axis_nodes_position_;
        x_axis_style = x_axis_style_;

        position();
        set_colors();

        if (ordered_gene_wise.length === 0) { return clear_fig(); }

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
        rect.transition(t_last).attr('y', d => y[d[y_attr]]).attr('height', rect_height-spacing).style('fill', d => colors(d[values]));
        ytre.filter(node => node.height === 0).transition(t_last).attr('y', d => y[d.data.id])
                                                                 .attr('x', y_axis_leaves_x)
                                                                 .style('text-anchor', (y_axis_leaves_position === 'before' ? 'end' : 'start'))
                                                                 .style('font-size', y_font_size)
                                                                 .attr('dy', y_font_size);
        ytre.filter(node => node.depth > 0 && node.height > 0).transition(t_last).attr('transform', d => 'translate('+(y_axis_nodes_x + d.y0)+','+d.x1+')rotate(-90)');
        ytre.filter(node => node.depth > 0 && node.height > 0).select('.ytre_box').transition(t_last).attr('width', d => d.x1 - d.x0).attr('height', d => d.y1 - d.y0).style('stroke-dasharray', d => (y_axis_nodes_position === 'before' ? pointing_right(d) : pointing_left(d)));
        ytre.filter(node => node.depth > 0 && node.height > 0).select('.ytre_clear').transition(t_last).attr('x', d => d.x1 - d.x0 - y_axis_nodes_x_width).attr('height', y_axis_nodes_x_width).attr('width', y_axis_nodes_x_width);
        ycat.transition(t_last).attr('y', x_axis_leaves_y).attr('x', d => y_category_x[d]).attr('dy', (x_axis_leaves_position === 'before' ? y_cat_font_size : 0))
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
        xtre.filter(node => node.depth > 0 && node.height > 0).select('.xtre_box').transition(t_last).attr('width', d => d.x1 - d.x0).attr('height', d => d.y1 - d.y0).style('stroke-dasharray', d => (x_axis_nodes_position === 'before' ? pointing_right(d) : pointing_left(d)));
        xtre.filter(node => node.depth > 0 && node.height > 0).select('.xtre_clear').transition(t_last).attr('x', d => d.x1 - d.x0 - x_axis_nodes_y_height).attr('height', x_axis_nodes_y_height).attr('width', x_axis_nodes_y_height);
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
            .on('click', (d) => (x_attr === 'gene_id' ? GeneCards(d.data.name) : null))
            .call(d3.drag().on('drag', drag_x).on('end', drag_x_end))
            .style('opacity', 0).transition(t_last).style('opacity', 1);

        xtre.enter()
            .filter(d => d.height > 0 && d.depth > 0)  // internal nodes
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
                .style('font-size', xtre_label_font_size)
                .attr('dy', xtre_label_font_size+spacing)
                .attr('dx', '0.2em')
            .select(function() { return this.parentNode; })
                .append('clipPath')
                .attr('id', d => 'clip-'+d.data.id)
                    .append('use')
                    .attr('xlink:href', d => '#rect-'+d.data.id)
                .select(function() { return this.parentNode; })
            .select(function() { return this.parentNode; })
                .append('rect')
                .attr('class', 'xtre_clear')
                .attr('id', d => 'clear-' + d.data.id)
                .attr('x', d => d.x1 - d.x0 - x_axis_nodes_y_height)
                .attr('height', d => x_axis_nodes_y_height)
                .attr('width', d => x_axis_nodes_y_height)
                .styles(clear_styles)
                .on('mouseover', function() { d3.select(this).style('opacity', 1)}).on('mouseout', function() { d3.select(this).style('opacity', 0)})
                .on('click', remove_node)
            .select(function() { return this.parentNode; })
            .style('opacity', 0).transition(t_last).style('opacity', 1);

        xcat.enter()
            .append('text')
            .attr('class', 'xcat')
            .attr('id', d => d)
            .text(d => d)
            .attr('x', y_axis_leaves_x)
            .attr('y', d => x_category_y[d])
            .style('font-size', x_cat_font_size)
            .attr('dy', x_cat_font_size)
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
            .on('click', (d) => (y_attr === 'gene_id' ? GeneCards(d.data.name) : null))
            .call(d3.drag().on('drag', drag_y).on('end', drag_y_end))
            .style('opacity', 0).transition(t_last).style('opacity', 1);

        ytre.enter()
            .filter(d => d.height > 0 && d.depth > 0)  // internal nodes
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
                .style('font-size', ytre_label_font_size)
                .attr('dy', ytre_label_font_size+spacing)
                .attr('dx', '0.2em')
            .select(function() { return this.parentNode; })
                .append('clipPath')
                .attr('id', d => 'clip-'+d.data.id)
                    .append('use')
                    .attr('xlink:href', d => '#rect-'+d.data.id)
                .select(function() { return this.parentNode; })
            .select(function() { return this.parentNode; })
                .append('rect')
                .attr('class', 'ytre_clear')
                .attr('id', d => 'clear-' + d.data.id)
                .attr('x', d => d.x1 - d.x0 - y_axis_nodes_x_width)
                .attr('height', d => y_axis_nodes_x_width)
                .attr('width', d => y_axis_nodes_x_width)
                .styles(clear_styles)
                .on('mouseover', function() { d3.select(this).style('opacity', 1)}).on('mouseout', function() { d3.select(this).style('opacity', 0)})
                .on('click', remove_node)
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
            .style('font-size', y_cat_font_size)
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
            .attr('fill', d => colors(d[values]))
            .style('opacity', 0).transition(t_last).style('opacity', 1);


        g.select('#rect_resizer').attr('cx', d3.max(Object.values(x))+rect_width).attr('cy', d3.max(Object.values(y))+rect_height);
        g.select('#xtre_resizer').attr('cx', d3.max(Object.values(x))+rect_width).attr('cy', x_axis_nodes_y + (x_axis_nodes_position === 'before' ? x_axis_nodes_y_height : x_axis_nodes_y_height*x_tree.height));
        g.select('#ytre_resizer').attr('cx', y_axis_nodes_x + (y_axis_nodes_position === 'before' ? y_axis_nodes_x_width : y_axis_nodes_x_width*y_tree.height)).attr('cy', d3.max(Object.values(y))+rect_height);


    }


    function set_colors() {
        values_domain = d3.extent(flatten(gene_wise), d => d[values]);

        if (color_style === 'interpolateTriplet') {
            colors = d3.scaleLinear().domain(values_domain.insert(1, 0)).range([negative_color, middle_color, positive_color]); }
        else {
            colors = d3.scaleSequential(d3[color_style]).domain(values_domain); }
    }

    function style({color_style_=color_style,
                    negative_color_=negative_color,
                    middle_color_=middle_color,
                    positive_color_=positive_color,
                    show_legends_=show_legends,
                    show_x_level_names_=show_x_level_names,
                    show_y_level_names_=show_y_level_names,
                    rotation_=rotation}={}) {

        color_style = color_style_;
        negative_color = negative_color_,
        middle_color = middle_color_,
        positive_color = positive_color_,
        show_legends = show_legends_;
        show_x_level_names = show_x_level_names_;
        show_y_level_names = show_y_level_names_;
        rotation = rotation_;

        // Colors
        set_colors();
        g.selectAll('.rect').style('fill', (d) => colors(d[values]));

        // Hide / Show & Rotate Labels
        x_axis_leaves_rotation = (x_axis_leaves_position === 'before') ? -rotation : rotation;

        g.selectAll('.xcat').attr('visibility', show_x_level_names ? 'visible' : 'hidden');
        g.selectAll('.ycat').attr('visibility', show_y_level_names ? 'visible' : 'hidden')
                            .attr('transform', d => 'rotate('+x_axis_leaves_rotation+','+y_category_x[d]+','+x_axis_leaves_y+')');

        g.selectAll('.xtre').filter(node => node.height === 0).attr('transform', d => 'rotate('+x_axis_leaves_rotation+','+d.x0+','+x_axis_leaves_y+')');

        // Legends
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
        ytre.filter(node => node.depth > 0 && node.height > 0).select('.ytre_label').style('font-size', ytre_label_font_size).attr('dy', ytre_label_font_size+spacing);
        ytre.filter(node => node.depth > 0 && node.height > 0).select('.ytre_clear').attr('x', d => d.x1 - d.x0 - y_axis_nodes_x_width).attr('height', y_axis_nodes_x_width).attr('width', y_axis_nodes_x_width);
        ycat.attr('y', x_axis_leaves_y)
            .attr('x', d => y_category_x[d])
            .attr('transform', d => 'rotate('+x_axis_leaves_rotation+','+y_category_x[d]+','+x_axis_leaves_y+')')
            .style('font-size', y_cat_font_size)
            .attr('dy', (x_axis_leaves_position === 'before' ? y_cat_font_size : 0));

        xtre.filter(node => node.height === 0).attr('x', d => d.x0)
                                              .attr('y', x_axis_leaves_y)
                                              .attr('transform', d => 'rotate('+x_axis_leaves_rotation+','+d.x0+','+x_axis_leaves_y+')')
                                              .style('font-size', x_font_size)
                                              .attr('dy', (x_axis_leaves_position === 'before' ? x_font_size : 0))
        xtre.filter(node => node.depth > 0 && node.height > 0).attr('transform', d => 'translate('+d.x0+','+(x_axis_nodes_y + d.y0)+')');
        xtre.filter(node => node.depth > 0 && node.height > 0).select('.xtre_box').attr('width', d => d.x1 - d.x0).attr('height', d => d.y1 - d.y0).style('stroke-dasharray', d => (y_axis_nodes_position === 'before' ? pointing_right(d) : pointing_left(d)));
        xtre.filter(node => node.depth > 0 && node.height > 0).select('.xtre_label').style('font-size', xtre_label_font_size).attr('dy', xtre_label_font_size+spacing);
        xtre.filter(node => node.depth > 0 && node.height > 0).select('.xtre_clear').attr('x', d => d.x1 - d.x0 - x_axis_nodes_y_height).attr('height', x_axis_nodes_y_height).attr('width', x_axis_nodes_y_height);
        xcat.attr('x', y_axis_leaves_x)
            .attr('y', d => x_category_y[d])
            .style('font-size', x_cat_font_size)
            .attr('dy', x_cat_font_size);

    }

    function clear_fig() {
        g.selectAll('.rect,.ytre,.xtre,.xcat,.ycat').transition(d3.transition().duration(500)).style('opacity', 0).remove();
    }

    ///////////////////////////////////////////////////////////////////////////
                          ///////      Drag      ///////
    ///////////////////////////////////////////////////////////////////////////

    function drag_node(d, hierarchy, xy, attr) {

        let initial_position = (node) => (xy === 'x' ? node.x0 : (node.height === 0 ? node.x0 : node.x1));

        index_of_dragging_node = d.data.order;
        current_position_of_dragging_node = (d.x ? d.x : initial_position(d)) + (xy === 'x' ? d3.event.dx : d3.event.dy);
        dragging_node_width = d.x1 - d.x0;

        set_nodes = d.parent.children;

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

        var nodes = g.selectAll('.'+xy+'tre').filter(node => node.data.id in updated_xy).each(node => {node.x = updated_xy[node.data.id]});

        g.selectAll('.rect').filter(rect => updated_xy[rect[attr]]).attr(xy, (rect) => updated_xy[rect[attr]]);
        nodes.filter(node => node.height === 0).attr(xy, node => node.x);
        nodes.filter(node => node.height > 0).attr('transform', (node) => {
            if (xy === 'x') { return 'translate('+node.x+','+(x_axis_nodes_y + node.y0)+')'; }
            if (xy === 'y') { return 'translate('+(y_axis_nodes_x+node.y0)+','+node.x+')rotate(-90)'; }
        });

        if (xy === 'x') {
            nodes.filter(node => node.height === 0).attr('transform', node => 'rotate('+x_axis_leaves_rotation+','+node.x+','+x_axis_leaves_y+')');
        }

    }

    function drag_node_end(d, hierarchy, xy, this_wise, other_wise) {

        set_nodes = d.parent.children;

        new_order = _.object(set_nodes.filter(node => node.data.id !== d.data.id)
                                       .map(node => [node.x, node.data.id])
                                       .concat([[d.x, d.data.id]])
                                       .sort((a, b) => a[0] - b[0])
                                       .map(([y, id], i) => [id, i]));

        hierarchy.each(node => { if (node.data.id && node.data.id in new_order) { node.data.order = new_order[node.data.id]} });
        hierarchy.each(node => { node.x = undefined; }); // do I even need this?

        old_index = hierarchy.leaves().map(leaf => leaf.data.id).indexOf(d.data.id);
        hierarchy.sort(function(a, b) { return a.depth - b.depth || a.data.order - b.data.order; });
        new_index = hierarchy.leaves().map(leaf => leaf.data.id).indexOf(d.data.id);

        this_wise.move(old_index, new_index);
        other_wise.forEach((other) => other.move(old_index, new_index));

        render();
        refresh_genes_cb();

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
            'x_axis_nodes_y_height_': Math.max(x_axis_nodes_y_height + (x_axis_nodes_position === 'before' ? -1 : 1) * (d3.event.dy / (x_tree.height-1)), 4),
            'rect_width_': Math.max((d3.event.x - x_tree.leaves().last().offset), x_tree.leaves().length*2) / x_tree.leaves().length,
        });
    }

    function drag_ytre_resizer(d) {
        d3.select(this).attr('cx', d3.event.x).attr('cy', d3.event.y);
        resize_fig({
            'y_axis_nodes_x_width_': Math.max(y_axis_nodes_x_width + (y_axis_nodes_position === 'before' ? -1 : 1) * (d3.event.dx / (y_tree.height-1)), 4),
            'rect_height_': Math.max((d3.event.y - y_tree.leaves().last().offset), y_tree.leaves().length*2) / y_tree.leaves().length,
        });
    }


    function remove_node(d) {

        if (d.ancestors().last().data.id === 'genes') {
            restart({'selected_gene_sets_': rendered_gene_sets().filter(gs => gs.gene_set_name !== d.data.name)});
            refresh_genes_cb();

        } else if (d.ancestors().last().data.id === 'metadata') {
            d.leaves().forEach(leaf => { delete samples_by_genes_matrix[leaf.data.name] });
            restart({'selected_gene_sets_': rendered_gene_sets()});
        }
    }

    function rendered_gene_sets() {
        if ('children' in genes) {
            return genes.children.map(node => {
                if (node.height === 0) { return {'gene_set_name':null, 'genes':[node.data.name]} }
                else { return {'gene_set_name': node.data.name, 'genes': node.children.map(gene => gene.data.name)} } });
        } else {
            return [];
        }
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
            if (t) {
                current_transform.x = clamp(-(x_tree.x1-100), w)(current_transform.x - d3.event.deltaY);
            } else {
                current_transform.y = clamp(-(y_tree.x1-100), h)(current_transform.y - d3.event.deltaY);
            }
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

        'rendered_gene_sets': rendered_gene_sets,

        transpose     : function() { t = !t; [rect_width, rect_height] = [rect_height, rect_width]; [y_axis_style, x_axis_style] = [x_axis_style, y_axis_style]; render(); },
        set_reordering: function(reordering_) { reordering = reordering_; if (reordering) { order(); render(); } },
    }

}





