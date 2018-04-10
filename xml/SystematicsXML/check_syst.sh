#!/bin/bash
echo -e "\n"
for filename in *.xml; do
	echo "Status and Checking of file" $filename ":"
	project.py --xml $filename --stage selection --status
	project.py --xml $filename --stage selection --checkana
	echo -e "\n"
done
