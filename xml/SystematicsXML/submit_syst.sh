#!/bin/bash
echo -e "\n"
for filename in *.xml; do
	echo "Cleaning and submitting file" $filename ":"
	project.py --xml $filename --stage selection --clean
	project.py --xml $filename --stage selection --submit
	echo -e "\n"
done
