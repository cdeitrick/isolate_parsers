from pathlib import Path


def main():
	class VCFColumns:
		sequence_id: str = 'seq id'
		position: str = 'position'
		alternate: str = 'alt'
		reference: str = 'ref'
		quality: str = 'quality'
		depth: str = 'readDepth'
		variant_type: str = 'variantType'

		def labels(self):
			for attribute in dir(self):
				if not attribute.startswith('_') and attribute != "labels":
					yield self.__getattribute__(attribute)
	VCFColumns = VCFColumns()
	from pprint import pprint
	for i in VCFColumns.labels():
		print(i)

if __name__ == "__main__":
	main()