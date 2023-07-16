# Cloning_verification_script
Cloning verification script that intakes barcodes and looks for matching sequence within barcodes based off target length.  To be used specifically for Crispr encoding plasmids/libraries that contain sgRNA sequence embedded in user-specified barcodes.  Detects MCS site and returns whether cloning was successful or not based off user inputted reference sheet.  Works with Sanger sequence output .ab1 files.

```
Adjustable settings found within main:
MCSregionSequence = 'GAAACACCGACTTGCAGGTGCTTAAGGGATCCAGTATACTGGATCGATTGATCACACCTGCGGAT'
reference_file = pd.read_excel('Reference.xlsx')
matching_column = 'Sequence_Construct'
first_barcode_sequence = 'ACCG'
final_barcode_sequence = 'GTTT'
target_length = 20
barcode_score_length = 8
approximate_matching = False
```
