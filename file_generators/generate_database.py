from pathlib import Path
import pandas
from pony.orm import Database, db_session, Optional, Required, PrimaryKey, Set
#annotation	description	gene	locusTag	mutationCategory	position	presentIn	presentInAllSamples	ref	seq id

class Mutation:
	id: PrimaryKey(auto = True)
	annotation = Optional(str)
	description = Optional(str)
	gene = Optional(str)
	locus_tag = Optional(str)
	category = Optional(str)
	seq_id = Required(str)
	position = Required(int)
	ref = Required(str)
	alt = Required(str)

class Sample:
	name: Required(str)
	samples: Set(Mutation)


def convert_table_to_database(comparison_table:Path, output_filename:Path):
	pass
if __name__ == "__main__":
	pass