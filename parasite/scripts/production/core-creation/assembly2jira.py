import argparse
import sys
import os
from ProductionUtils import print_info,exit_with_error,api_request
import re
import json
import getpass
from jira import JIRA

# NCBI Configuration
ncbi_accession_pattern = r'^(GCA|GCF)_\d+(\.\d+)?$'
ncbi_api_url = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{0}/dataset_report"

# JIRA Configuration
jira_project_key = "PARASITE"
jira_subtask_issuetype_id = "5"
jira_server = 'https://www.ebi.ac.uk/panda/jira/'
jira_parent_key = os.environ["PARASITE_JIRA_PREPROCESSING_KEY"]

# Configuration
helpdesk_url_for_id = "https://helpdesk.ebi.ac.uk/Ticket/Display.html?id={0}"
parasite_version = os.environ["PARASITE_VERSION"]

# Either species or ncbi_assembly_acc must be provided by the user
class SpeciesAssemblyAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if not (namespace.species or namespace.ncbi_assembly_acc):
            parser.error("Either --species or --ncbi_assembly_acc must be provided.")
        setattr(namespace, self.dest, values)

# Format the help message
class CustomFormatter(argparse.HelpFormatter):
    def _format_usage(self, usage, actions, groups, prefix):
        examples = [
            'Example Usage:\n python ${WORM_CODE}/parasite/scripts/production/core-creation/assembly2jira.py \
--ncbi_assembly_acc GCA_030248215.1 --type "Assembly Update" --name "Jason Tsai" --email "ijtsai@gate.sinica.edu.tw" --helpdesk 655052 \
--notes "Data deposited [here|https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/datasets/PRJNA953805/enoplolaimus_lenunculus/PRJNA953805]"\n'
        ]
        usage_str = super()._format_usage(usage, actions, groups, prefix)
        return f"{usage_str}\n\n{' '.join(examples)}\n"

def get_strain(assembly_summary):
    strain_values = []
    try:
        strain_values.append(assembly_summary['organism']['infraspecific_names']['strain'])
    except: pass
    try:
        for sinfo in assembly_summary['assembly_info']['biosample']['attributes']:
            if 'name' in sinfo.keys() and 'value' in sinfo.keys():
                if sinfo['name']=='strain' and sinfo['value']!='not applicable':
                    strain_values.append(sinfo['value'])
                if sinfo['name']=='isolate' and sinfo['value']!='not applicable':
                    strain_values.append(sinfo['value'])
    except:
        pass
    if len(strain_values)>0:
        return ",".join(list(set(strain_values)))
    else:
        return ""
    
def main():
    parser = argparse.ArgumentParser(description="Script that will accept an INSDC assembly accession and will \
                                     create a JIRA ticket for the core creation process. The ticket will be \
                                     created as a subtask of the PARASITE_JIRA_PREPROCESSING_KEY ticket.",
                                     formatter_class=CustomFormatter)
    parser.add_argument('--species', help='Species name', default=None)
    parser.add_argument('--bioproject', help='BioProject ID (Optional)')
    parser.add_argument('--insdc_assembly_name', help='INSDC Assembly name (Optional)')
    parser.add_argument('--ncbi_assembly_acc', help='NCBI Assembly accession (Optional)')
    parser.add_argument('--strain_biosample', help='Strain or BioSample (Optional)')
    parser.add_argument('--publication', help='Publication (Optional)')
    parser.add_argument('--type', required=True, help='Type (Required)')
    parser.add_argument('--status', default='READY', help='Status (Optional, Default=READY)')
    parser.add_argument('--notes', help='Notes (Optional)')
    parser.add_argument('--name', help='Name (Optional)')
    parser.add_argument('--email', help='Email (Optional)')
    parser.add_argument('--helpdesk', help='Helpdesk (Optional)', action=SpeciesAssemblyAction)


    args = parser.parse_args()
    
    # Process the command-line arguments here and perform the desired actions
    # For example, you can access the options using args.species, args.bioproject, etc.
    # Print or use the values as per your requirement
    
    print_info("Arguments:")
    if args.species: print_info("Species: " + args.species)
    if args.bioproject: print_info("BioProject: " + args.bioproject)
    if args.insdc_assembly_name: print_info("INSDC Assembly Name: " + args.insdc_assembly_name)
    if args.ncbi_assembly_acc: print_info("NCBI Assembly Acc: " + args.ncbi_assembly_acc)
    if args.strain_biosample:print_info("Strain/BioSample: " + args.strain_biosample)
    if args.publication: print_info("Publication: " + args.publication)
    if args.type: print_info("Type: " + args.type)
    if args.status: print_info("Status: " + args.status)
    if args.notes: print_info("Notes: " + args.notes)
    if args.name: print_info("Name: " + args.name)
    if args.email: print_info("Email: " + args.email)
    if args.helpdesk: print_info("Helpdesk: " + args.helpdesk)

    if args.ncbi_assembly_acc:
        if re.match(ncbi_accession_pattern, args.ncbi_assembly_acc) is None:
            exit_with_error(f"Invalid NCBI Assembly Accession: {args.ncbi_assembly_acc}")
        ncbi_assembly_acc = args.ncbi_assembly_acc
        assembly_summary_response = json.loads(api_request(ncbi_api_url.format(ncbi_assembly_acc)))
    if assembly_summary_response['total_count']>1 or assembly_summary_response['total_count']<1:
        print(assembly_summary_response)
        print(f"NCBI API Request for {args.ncbi_assembly_acc} returned multiple or no reports.")

    assembly_summary = assembly_summary_response["reports"][0]

    if args.species:
        species = args.species
    else:
        species = assembly_summary['organism']['organism_name']
    
    if args.bioproject:
        bioproject = args.bioproject
    else:
        bioproject = assembly_summary['assembly_info']['bioproject_accession']
        
    if args.insdc_assembly_name:
        insdc_assembly_name = args.insdc_assembly_name
    else:    
        insdc_assembly_name = assembly_summary['assembly_info']['assembly_name']
    
    if args.strain_biosample:
        biosample = args.strain_biosample
    else:
        biosample = assembly_summary['assembly_info']['biosample']['accession']
    
    if args.strain_biosample:
        strain_biosample = args.strain_biosample
    else:
        strain_biosample = f"{get_strain(assembly_summary)}/{biosample}"

    genome = species.lower().replace(" ","_") + "_" + bioproject.lower()
    helpdesk_url = helpdesk_url_for_id.format(args.helpdesk)

    notes = "Submitted by: " + assembly_summary['assembly_info']['submitter'] + ","
    notes += ("Submitted by: " + f"[{args.name}|mailto:{args.email}], " if args.email and args.name else "")
    notes += ("Submitted by: " + f"{args.name}, " if args.name and not args.email else "")
    notes += ("Submitted by: " + f"{args.email}, " if args.email and not args.name else "")
    notes += ("Submitted by: " + f"[{args.email}|mailto:{args.email}], " if args.email and not args.name else "")
    notes += (f"\nHelpdesk ticket: [{args.helpdesk}|{helpdesk_url}]." if args.helpdesk else "")
    notes += (args.notes if args.notes else "")

    username = getpass.getuser()
    password = getpass.getpass(f"Enter your Jira password for {username}: ")
    jira = JIRA(server=jira_server,basic_auth=(username, password))

    #create subtask for parent_issue
    # create the subtask
    subtask = jira.create_issue(
        project = jira_project_key,
        summary = f"WBPS{parasite_version} - {genome}",
        description = f"|{species}|{bioproject}|{insdc_assembly_name}|{ncbi_assembly_acc}|{strain_biosample}|" + \
                                            (f"{args.publication}|" if args.publication else "") + \
                                            f"{args.type}|{args.status}|{notes}|",
        issuetype={'name': 'Sub-task'},
        parent={'key': jira_parent_key}
    )

    subtask_browse_url = jira_server + f"browse/{subtask.key}"
    subtask_browse_link = f"[{subtask.key}|{subtask_browse_url}]"

    # modify subtask description to include link to subtask itself
    subtask.update(description=f"|{species}|{subtask_browse_link}|{bioproject}|{insdc_assembly_name}|{ncbi_assembly_acc}|{strain_biosample}|" + \
                                            (f"{args.publication}|" if args.publication else "") + \
                                            f"{args.type}|{args.status}|{notes}|")


    # add the link to the parent issue
    jira.add_comment(jira_parent_key, f"Subtask created: {subtask_browse_link}")

    print_info(f"Subtask created: {subtask}")
    print_info(f"URL: {jira_server}/browse/{subtask.key}")
    print_info(f"\n\nDone")


    
if __name__ == "__main__":
    main()
