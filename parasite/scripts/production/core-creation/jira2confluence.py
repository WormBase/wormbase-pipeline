from jira import JIRA
import getpass
import argparse
import os
import re

parser = argparse.ArgumentParser(description="Script that will accept a PARASITE " +\
    "genome preprocessing JIRA ticket key and will output a TSV will all the subtask " +\
    "descriptions. The output TSV file can then be pasted into Confluence Release tables " +\
    "here: https://www.ebi.ac.uk/seqdb/confluence/display/WORMBASE/ParaSite+genome+list")   
parser.add_argument('--jira_mother_key', help='Mother pre-processing JIRA key ID. Example: PARASITE-529', default=os.environ["PARASITE_JIRA_PREPROCESSING_KEY"])

args = parser.parse_args()

def is_jira_link(input_string):
    # Find the indices of "[" and "]"
    open_bracket_index = input_string.find("[")
    close_bracket_index = input_string.find("]")
    # Check if "[" and "]" exist and "[" comes before "]"
    if open_bracket_index != -1 and close_bracket_index != -1 and open_bracket_index < close_bracket_index:
        # Get the substring between "[" and "]"
        substring_between_brackets = input_string[open_bracket_index + 1 : close_bracket_index]
        # Check if "|" is present in the substring
        if "|" in substring_between_brackets:
            return True
    return False

PARASITE_SCRATCH = os.environ["PARASITE_SCRATCH"]

# JIRA Configuration
jira_project_key = "PARASITE"
jira_subtask_issuetype_id = "5"
jira_server = 'https://www.ebi.ac.uk/panda/jira/'
jira_parent_key = args.jira_mother_key
clean_description_pattern = r'\|\*Species\*\|.+\|\*Notes\*\|'
link_pattern= r'(\[[^\]]+\|[^\]]+\])'

# Output Configuration
output_dir = os.path.join(PARASITE_SCRATCH,"jira")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
output_file = os.path.join(output_dir,f'{jira_parent_key}_descriptions.tsv')

# Initialize the JIRA client
# username = getpass.getuser()
# password = getpass.getpass(f"Enter your Jira password for {username}: ")
username = getpass.getuser()
password = getpass.getpass(f"Enter your Jira password for {username}: ")
jira = JIRA(server=jira_server,basic_auth=(username, password))

# Parent issue key
parent_issue_key = 'PARASITE-529'

# Get the parent issue
parent_issue = jira.issue(parent_issue_key)

# Get the subtasks of the parent issue
subtasks = jira.search_issues(f'parent = {parent_issue.key}')

# Gather descriptions of all subtasks
all_descriptions = []
for subtask in subtasks:
    description = subtask.fields.description
    clean_description1 = re.sub(clean_description_pattern, '', description).strip()
    clean_description2 = re.sub("{color(:[^}]+)?}","",clean_description1)
    split_links = re.split(link_pattern,clean_description2)
    split_links2 = [re.sub(r'\s+', ' ', x) for x in split_links]
    clean_description3 = "".join([x if is_jira_link(x) else x.replace("|","\t") for x in split_links2])
    clean_description4 = re.sub(r'(\s+)?\\t(\s+)?', '\t', clean_description3)
    clean_description5 = clean_description4.replace(u'\xa0', u' ')
    clean_description6 = re.sub(r' +', ' ', clean_description5)
    clean_description7 = clean_description6.strip("\t").strip().replace("\n"," ")
    all_descriptions.append(clean_description7)

# Write the descriptions to a TSV file
with open(output_file, 'w') as f:
    for description in all_descriptions:
        f.write(description + '\n')
