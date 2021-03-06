#! /bin/sh

set -e

usage () {
    echo "atom-autocomplete-clang-gen [project path]
generates a .clang_autocomplete file including all paths containing c/c++
headers. additional options can be provided in a .clang_autocomplete.addon file
present in the project dir. Special values are allowed in this addon:
- {PROJECT_PATH} : Full path of the project path
"

}

if [ $# -ne 1 ] ; then
   usage
   exit -1
fi

project_path=$1
config_file=${project_path}/.clang_complete
config_addon_file=${project_path}/.clang_complete.addon

# Add project dirs containing include files
include_path_list=`find -L $project_path -name "*.h" | xargs -I{} dirname {} | uniq`
rm -f $config_file
for path in $include_path_list
do
    echo "-I$path" >> $config_file
done

# Add addon config
escaped_project_path=`echo $project_path | sed "s,/,\\\\/,g"`
if [ -f "$config_addon_file" ]
then
    cat $config_addon_file | sed "s,{PROJECT_PATH},$escaped_project_path,g"  >> $config_file
fi


echo -e "\033[32m[Done]\033[0m"
