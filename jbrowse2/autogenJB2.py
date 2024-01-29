import os

from jbrowse2 import jbrowseConnector




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="", epilog="")
    parser.add_argument("--yaml", help="Track Configuration")
    parser.add_argument("--outdir", help="Output directory", default="out")
    parser.add_argument("--version", "-V", action="version", version="%(prog)s 0.0.1")
    args = parser.parse_args()
    if os.path.exists(args.outdir):
        root = args.outdir
        dirList =  os.scandir(root)
        subdirs = [f.path for f in dirList if f.is_dir()]
        genome_paths = [f.path for f in dirList if f.name.startswith('genome') and f.is_file()]
        if len(genome_paths) > 0:
            genome_fnames = [os.path.basename(x).split('_')[2:] for x in genome_paths]  # expect genome_1_genomename.fasta etc
            jc = JbrowseConnector(
                outdir=args.outdir,
                genomes=[
                    {
                        "path":  x,
                        "meta": {"name" : genome_fnames[i], },
                    }
                    for i,x in enumerate(genome_paths)
                ],
            )
            jc.process_genomes()
            # .add_default_view() replace from https://github.com/abretaud/tools-iuc/blob/jbrowse2/tools/jbrowse2/jbrowse2.py
            default_session_data = {
                "visibility": {
                    "default_on": [],
                    "default_off": [],
                },
                "style": {},
                "style_labels": {},
            }

            track_paths = [x for x in genome_paths if not x.startswith('genome') and  x.is_file()]
            for i, track in enumerate(track_paths):
                track_conf = {}
                track_conf['format']= os.path.basename(track).split('_')[0]
                track_conf["name"] = os.path.basename(track).split('_')[2:]   # expect genome_1_genomename.fasta etc
                fext = os.path.splitext(os.path.basename(track)).replace('.','')
                track_conf["label"] = "%s_%i" % (os.path.basename(track), i)
                track_conf["trackfiles"] = []
                keys = jc.process_annotations(track_conf)

                if keys:
                    for key in keys:
                        default_session_data["visibility"][
                            track.attrib.get("visibility", "default_off")
                        ].append(key)
                        if track_conf.get("style", None):
                            default_session_data["style"][key] = track_conf[
                                "style"
                            ]  # TODO do we need this anymore?
                        if track_conf.get("style_lables", None):
                            default_session_data["style_labels"][key] = track_conf.get(
                                "style_labels", None
                            )
            # default_session_data["defaultLocation"] = root.find(
                # "metadata/general/defaultLocation"
            # ).text
            # default_session_data["session_name"] = root.find(
                # "metadata/general/session_name"
            # ).text
            # general_data = {
                # "analytics": root.find("metadata/general/analytics").text,
                # "primary_color": root.find("metadata/general/primary_color").text,
                # "secondary_color": root.find("metadata/general/secondary_color").text,
                # "tertiary_color": root.find("metadata/general/tertiary_color").text,
                # "quaternary_color": root.find("metadata/general/quaternary_color").text,
                # "font_size": root.find("metadata/general/font_size").text,
            # }
            # jc.add_general_configuration(general_data)
            trackconf = jc.config_json.get("tracks", None)
            if trackconf:
                jc.config_json["tracks"].update(jc.tracksToAdd)
            else:
                jc.config_json["tracks"] = jc.tracksToAdd
            jc.write_config()
            jc.add_default_session(default_session_data)
            # jc.text_index() not sure what broke here.
