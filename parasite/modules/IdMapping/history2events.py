from ProductionMysql import *
from ProductionUtils import *
import csv
# parser = OptionParser(usage='usage: %prog [options] arguments')
# parser.add_option("-c", "--core_db", dest="COREDB",
#                   help="Required: core db name pattern")
# parser.add_option("-f", "--from", dest="FROM",
#                   help="Required: from assembly (old_assembly) name")
# parser.add_option("-f", "--from", dest="TO",
#                   help="Required: to assembly (new_assembly) name")
#
# (options, args) = parser.parse_args()
#
# COREDB = options.COREDB.strip()
# FROM = options.FROM.strip()
# TO = options.TO.strip()

CORE_DB = "mansoni"
FROM = "ASM23792v2"
TO = "SM_V9"

core = Core(staging.host, "mansoni")
mapping_sessions = {x[0]:{"old_assembly":x[1], "new_assembly":x[2]}
                    for x in core.connect().execute("select mapping_session_id,old_assembly,new_assembly from mapping_session;")}

ste = {x:{"old_genes":list(set([x[0] for x in core.connect().execute("select old_stable_id from stable_id_event where mapping_session_id="+str(x))])),
          "new_genes":list(set([x[0] for x in core.connect().execute("select new_stable_id from stable_id_event where mapping_session_id="+str(x))]))}
       for x in mapping_sessions}

all_genes = [x[0] for x in core.connect().execute("SELECT old_stable_id FROM stable_id_event UNION (SELECT new_stable_id FROM stable_id_event)")]
latest_genes = [x[0] for x in core.connect().execute("SELECT stable_id from gene")]

def old_gene_per_session(gene_id, session):
    fres = [list(x) for x in core.connect().execute("select * from stable_id_event where mapping_session_id="+str(session)+
                                       " and old_stable_id="+"'"+gene_id+"'")]
    return(fres)

def new_gene_per_old_gene_per_session(gene_id, session):
    fres = [x for x in core.connect().execute("select new_stable_id from stable_id_event where mapping_session_id="+str(session)+
                                       " and old_stable_id="+"'"+gene_id+"'")]
    return(fres)

def new_gene_per_session(gene_id, session):
    fres = [x for x in core.connect().execute("select * from stable_id_event where mapping_session_id="+str(session)+
                                       " and new_stable_id="+"'"+gene_id+"'")]
    return(fres)

def old_gene_per_session_if_new_null(gene_id, session):
    fres = [x for x in core.connect().execute("select * from stable_id_event where mapping_session_id="+str(session)+
                                       " and old_stable_id="+"'"+gene_id+"'"+" and new_stable_id is null")]
    return(fres)

def new_gene_per_session_if_old_null(gene_id, session):
    fres = [x for x in core.connect().execute("select * from stable_id_event where mapping_session_id="+str(session)+
                                       " and new_stable_id="+"'"+gene_id+"'"+" and old_stable_id is null")]
    return(fres)

final_dict = {}
for gene in all_genes:
    if gene is not None:
        if gene in ste[list(mapping_sessions.keys())[0]]["old_genes"]:
            if gene in latest_genes:
                print(gene+" preserved")
                final_dict[gene] = {"SM_V5_ID":gene, "SM_V9_ID":gene, "status":"preserved"}
            elif gene not in latest_genes:
                #print(gene)
                if max([len(old_gene_per_session(gene, mid)) for mid in list(mapping_sessions.keys())])>1:
                    print(gene + " 1-to-many")
                    final_dict[gene] = {"SM_V5_ID": gene, "SM_V9_ID": "None", "status": "1-to-many"}
                elif max([len(old_gene_per_session_if_new_null(gene, mid)) for mid in list(mapping_sessions.keys())])==1:
                    print(gene + " deleted")
                    final_dict[gene] = {"SM_V5_ID": gene, "SM_V9_ID": "None", "status": "deleted"}
                elif max([len(old_gene_per_session(gene, mid)) for mid in list(mapping_sessions.keys())])==1:
                    new_genes = flatten([x[0] for x in [new_gene_per_old_gene_per_session(gene, mid) for mid in list(mapping_sessions.keys())] if len(x)!=0])
                    if len(new_genes) == 1:
                        if max([len(new_gene_per_session(new_genes[0], mid)) for mid in list(mapping_sessions.keys())]) > 1:
                            print(gene + " merged to " + new_genes[0])
                            final_dict[gene] = {"SM_V5_ID": gene, "SM_V9_ID": new_genes[0], "status": "merged_to:"+new_genes[0]}
                        else:
                            print(gene + " renamed to " + new_genes[0])
                            final_dict[gene] = {"SM_V5_ID": gene, "SM_V9_ID": new_genes[0], "status": "renamed_to:"+new_genes[0]}
                    else:
                        print("ERROR: how many times has "+gene+" been merged/renamed? Investigate.")
                        final_dict[gene] = {"SM_V5_ID": "None", "SM_V9_ID": "None", "status": "ERROR"}
                else:
                    print("ERROR: Don't know what to do with: "+gene)
                    final_dict[gene] = {"SM_V5_ID": "None", "SM_V9_ID": "None", "status": "ERROR"}
        elif gene not in ste[list(mapping_sessions.keys())[0]]["old_genes"]:
            if gene in latest_genes:
                if max([len(new_gene_per_session_if_old_null(gene, mid)) for mid in list(mapping_sessions.keys())[1:]]) == 1:
                    print(gene + " new gene")
                    final_dict[gene] = {"SM_V5_ID": "None", "SM_V9_ID": gene, "status": "new"}
                else:
                    max_new_genes_per_session = max([len(new_gene_per_session(gene, mid)) for mid in list(mapping_sessions.keys())])
                    if  max_new_genes_per_session == 1:
                        print(gene + " new gene (production of a 1-to-many situation)")
                        final_dict[gene] = {"SM_V5_ID": "None", "SM_V9_ID": gene, "status": "new_from_1tomany"}
                    elif max_new_genes_per_session > 1:
                        print(gene + " new gene (production of a merge)")
                        final_dict[gene] = {"SM_V5_ID": "None", "SM_V9_ID": gene, "status": "new_from_merge"}
                    else:
                        print("ERROR: Don't know what to do with: "+gene)
                        final_dict[gene] = {"SM_V5_ID": "None", "SM_V9_ID": "None", "status": "ERROR"}
            else:
                print(gene + " deleted")
                final_dict[gene] = {"SM_V5_ID": "None", "SM_V9_ID": "None", "status": "deleted"}
        else:
            print("ERROR: Don't know what to do with: " + gene)
            final_dict[gene] = {"SM_V5_ID": "None", "SM_V9_ID": "None", "status": "ERROR"}

with open('/homes/digri/smansoni_SM5_vs_SM9.tsv', 'w') as f:
    f.write("gene_ID\tSM_V5_gene_ID\tSM_V9_gene_ID\tstatus\n")
    for key in final_dict.keys():
        f.write("%s\t%s\t%s\t%s\n" % (key, final_dict[key]["SM_V5_ID"], final_dict[key]["SM_V9_ID"], final_dict[key]["status"]))


