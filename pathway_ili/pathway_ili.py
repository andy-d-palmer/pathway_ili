from flask import Flask, render_template, request, send_file

from .utils import generate_ili_csv, get_map_filename

app = Flask(__name__)

@app.route('/')
def index():
    return render_template('index.html')


@app.route('/pathway/<map_id>/ili/')
def get_ili_map(map_id):
    group_1_name = request.args.get('group_1_name', 'group1', type=str)
    group_2_name = request.args.get('group_2_name', 'group2', type=str)
    group_3_name = request.args.get('group_3_name', 'group3', type=str)
    group_1_ds_names = request.args.get('group_1_ds_names', '', type=str)
    group_2_ds_names = request.args.get('group_2_ds_names', '', type=str)
    group_3_ds_names = request.args.get('group_3_ds_names', '', type=str)
    print (group_1_name, group_1_ds_names)
    print (group_2_name, group_2_ds_names)
    print (group_3_name, group_3_ds_names)

    groups = {}
    for gn, g in zip([group_1_name, group_2_name, group_3_name], [group_1_ds_names, group_2_ds_names, group_3_ds_names]):
        if not g == '':
            groups[gn] = g.split("|")
    ili_csv = generate_ili_csv(groups, map_id)

    return send_file(ili_csv,
                     mimetype='text/csv',
                     attachment_filename="ili_{}.csv".format(map_id),
                     as_attachment=True
                     )

@app.route('/pathway/<map_id>/im/')
def get_map(map_id):
    filename = get_map_filename(map_id)
    return send_file(filename,
                     mimetype='image/png',
                     attachment_filename="{}.png".format(map_id),
                    as_attachment=True)

if __name__ == '__main__':
    app.run()