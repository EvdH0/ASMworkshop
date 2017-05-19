## ASM Microbe '17 workshop
Room 365

# Timetable
- Morten Sommer:  Introduction to antibiotic resistome
- Lejla Imamovic:  Resistome profiling 
(from the lab to databases, overview of the data and analysis to be used during the workshop)
- Mostafa Ellabaan (with Eric’s & Lejla’s help):
Task 1.        Introduce basic command line and tools how to get the data on/off server
      Output 1:   Ability to use basic commands needed for the workshop (data, database locations etc)

Task 2.        Assign ORFs in GeneMark for selected dataset ((MinION last run)
      Output 2:    Protein (pro) and nucleotide (nt) sequences in fasta format
      Question:    How many ORFs are predicted?
      Task 3.        Format databases (e.g. CARD nt and pro or Resfinder nucletide)
      Output 3:    Database for further use, knowledge on importance of databases
                           Also mention  that one can curate their own database  
                           Task 4A.     Annotate insert nt sequences with Resfinder online server (promote DTU)
Task 4B.      Annotate ORFs nt sequences against Resfinder database (blastn)
      Output 4:    Manually curated horizontally transferred list of genes from inserts/contigs  and ORFs. 
      Questions:  How many inserts have antibiotic resistance gene?
                            Is there a difference between DTU web server (contig) (4A) and annotated ORFs (4B)?
Task 5.       Annotate protein sequences against CARD protein database
      Output 5:   List of antibiotic resistance genes from ORFs. 
      Questions: How many inserts have antibiotic resistance gene? 
                           Input from participant of interesting hits?

Task 6.       Annotate protein sequences using hmm model in Pfam 
      Output 6:   List of genes from ORFs. 
      Questions: What can we see from Pfam domain list?
If there is time left, annotate ORFs against optimal databases using plasmid database (to look for mobility). 
We can use megablast to map the nt ORFs to plasmid database.

- Eric van der Helm:  Rapid resistome profiling using Nanopore sequencing
- Morten Sommer:  Antibiotic resistance gene exchange networks

# Workshop

## Login in to the amazon cloud

## Prepare your directories

## Finding Open Reading Frames (ORFs)

## The first database, Resfinder

## CARD

## Combine data






You can use the [editor on GitHub](https://github.com/EvdH0/ASM/edit/master/README.md) to maintain and preview the content for your website in Markdown files.

Whenever you commit to this repository, GitHub Pages will run [Jekyll](https://jekyllrb.com/) to rebuild the pages in your site, from the content in your Markdown files.

### Markdown

Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/EvdH0/ASM/settings). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://help.github.com/categories/github-pages-basics/) or [contact support](https://github.com/contact) and we’ll help you sort it out.
