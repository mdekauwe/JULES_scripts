#!/bin/bash

from="trunk@20512"
to="vn6.1_trunk"
svn co https://code.metoffice.gov.uk/svn/jules/main/$from $to

cd $to

cp ../../JULES_scripts/compile_jules_on_mac_wrapper.sh .
