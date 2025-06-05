import argparse
import json
import logging
import sys
from typing import List, Tuple

import pysam

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class VCFFilterError(Exception):
    pass

class JSONCriteriaParser:
    VALID_OPERATORS = ['>=', '<=', '>', '<', '==', '!=']

    def __init__(self, criteria_file: str):
        self.criteria = self._load_criteria(criteria_file)
        self.parsed_criteria = self._parse_criteria()
        self.multiallelic_strategy = self.criteria.get('_multiallelic_strategy', 'any')

    def _load_criteria(self, criteria_file):
        try:
            with open(criteria_file) as f:
                criteria = json.load(f)
            if not isinstance(criteria, dict):
                raise VCFFilterError("Criteria file must be a JSON object")
            return criteria
        except Exception as e:
            raise VCFFilterError(f"Failed to load criteria: {e}")

    def _parse_criteria(self) -> List[Tuple[str, str, str]]:
        parsed = []
        for field, condition in self.criteria.items():
            if field.startswith('_'):
                continue
            for op in sorted(self.VALID_OPERATORS, key=len, reverse=True):
                if op in condition:
                    lhs, rhs = condition.split(op, 1)
                    parsed.append((field, op, rhs.strip()))
                    break
            else:
                raise VCFFilterError(f"Invalid condition for field '{field}': {condition}")
        return parsed

class VCFProcessor:
    def __init__(self, input_file, output_file, criteria_parser):
        self.input_file = input_file
        self.output_file = output_file
        self.criteria_parser = criteria_parser
        self.processed = 0
        self.passed = 0

    def process(self):
        try:
            in_vcf = pysam.VariantFile(self.input_file)
            out_vcf = pysam.VariantFile(self.output_file, 'w', header=in_vcf.header)

            for rec in in_vcf:
                self.processed += 1
                if self._passes_criteria(rec):
                    rec.filter.clear()
                    rec.filter.add("PASS")
                    self.passed += 1
                out_vcf.write(rec)

            logger.info(f"Processed {self.processed} records, {self.passed} passed.")
        except Exception as e:
            raise VCFFilterError(f"Error processing VCF: {e}")

    def _passes_criteria(self, rec):
        for field, op, expected in self.criteria_parser.parsed_criteria:
            val = rec.info.get(field, None)
            if val is None:
                return False

            if isinstance(val, tuple):
                val = list(val)
            elif not isinstance(val, list):
                val = [val]

            strategy = self.criteria_parser.multiallelic_strategy
            if strategy == 'any':
                return any(self._compare(v, op, expected) for v in val)
            elif strategy == 'all':
                return all(self._compare(v, op, expected) for v in val)
            elif strategy == 'first':
                return self._compare(val[0], op, expected)
            else:
                return any(self._compare(v, op, expected) for v in val)
        return True

    def _compare(self, actual, op, expected):
        try:
            actual = float(actual)
            expected = float(expected)
        except ValueError:
            pass

        if op == '==': return actual == expected
        if op == '!=': return actual != expected
        if op == '>=': return actual >= expected
        if op == '<=': return actual <= expected
        if op == '>': return actual > expected
        if op == '<': return actual < expected
        return False

def main():
    parser = argparse.ArgumentParser(description="Filter VCF using criteria from JSON")
    parser.add_argument("-i", "--input", required=True, help="Input VCF")
    parser.add_argument("-o", "--output", required=True, help="Output VCF")
    parser.add_argument("-c", "--criteria", required=True, help="Criteria JSON")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        criteria = JSONCriteriaParser(args.criteria)
        processor = VCFProcessor(args.input, args.output, criteria)
        processor.process()
    except VCFFilterError as e:
        logger.error(e)
        sys.exit(1)

if __name__ == '__main__':
    main()
