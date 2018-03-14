#!/bin/bash

export PATH="/home/users/kjyi/tools/MToolBox/MToolBox:$PATH"

cat <<EOF > .05a_config.sh
#!/bin/bash
input_path=/home/users/kjyi/Projects/chromophobe/.fastq/
output_name=/home/users/kjyi/Projects/chromophobe/mtoolbox
input_type=fastq
ref=RSRS
UseMarkDuplicates=false
UseIndelRealigner=false
MitoExtraction=false
EOF

MToolBox.sh -i /home/users/kjyi/Projects/chromophobe/.05a_config.sh -m "-t 20" -a "-z 0.1"
