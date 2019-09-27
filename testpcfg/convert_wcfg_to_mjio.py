#convert_trees_to_pcfg.py

import utility
import wcfg
import argparse


parser = argparse.ArgumentParser(description='Convert Grammar from WCFG format to MJIO')
parser.add_argument('input', type=str, help='filename of input grammar')
parser.add_argument('output', type=str, help='filename of output grammar')

args = parser.parse_args()

mywcfg = wcfg.load_wcfg_from_file(args.input)

# Convert to pcfg 

mywcfg.renormalise()
mywcfg.store_mjio(args.output)

