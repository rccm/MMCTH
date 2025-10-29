# demo_args.py
import argparse
p = argparse.ArgumentParser()
p.add_argument("--name", required=True)
p.add_argument("--age", type=int, required=True)
args = p.parse_args()
print(f"Hello {args.name}, next year you'll be {args.age + 1}.")