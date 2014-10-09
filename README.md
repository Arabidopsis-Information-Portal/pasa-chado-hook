pasa-chado-hook
===============

PASA Data Adapter to interact with a Chado Schema database

* Annot_retriever extracts gene structures from Chado, generates an intermediate `Gene_obj` index before loading into PASA database
    * STATUS: Complete
* Annot_updater processes gene structure updates present in the PASA database (new, merge, split, alt_splice) and writes back to Chado
    * STATUS: Under development

Required Perl modules:
* Carp
* DBI
* DBD:Pg
* FindBin
* URI::Escape

Required PASA-specific Perl modules:
* Gene_obj
* Gene_obj_indexer

Reference material:

* Writing PASA data adapters: [http://pasa.sourceforge.net/PASA_data_adapters.html]()
* Chado `Sequence` module: [http://gmod.org/wiki/Chado_Sequence_Module]()

