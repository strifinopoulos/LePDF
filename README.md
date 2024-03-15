
# LePDF

## Standard Model PDFs of Leptons

Authors: *Francesco Garosi, David Marzocca, Sokratis Trifinopoulos*

Reference: [https://arxiv.org/abs/2303.16964]

------------------

### Description

This repository contains files for the parton distribution functions (PDFs) of electrons, muons, as well as the corresponding antiparticles. These have been obtained by solving numerically the DGLAP equations starting from the scale of the lepton mass up to several tens of TeV. All Standard Model interactions are included, starting from the corresponding relevant mass scale.
For more details on the derivation see the reference linked above and for instructions on how to use LePDF or on the source code implementation see *LePDF_manual.pdf*.


### PDF sets

The repository contains 8 different PDF sets, each named as *LePDF_(particle)_(scheme)_0000.dat*, with the corresponding .info file *LePDF_(particle)_(scheme).info*.

The field *(particle)* can be: *e* (electron), *eb* (positron), *mu* (muon), *mub* (anti-muon).

The field *(scheme)* can be either *5FS* (5-flavour scheme) or *6FS* (6-flavour scheme), in which case the top quark is not (or is) included in the evolution.

The PDF sets are exported in a format very similar to *LHAPDF6* [https://arxiv.org/abs/1412.7420], which is however modified by adding, for each particle, a further specification for the different helicity state, since PDFs become polarised above the electroweak scale. For more details on the format see Appendix F of [https://arxiv.org/abs/2303.16964] or *LePDF_manual.pdf*.

### Examples

With the release we also include a Mathematica notebook, *LePDF_examples.nb*, which includes examples of how these PDF sets can be loaded and used to produce plots, parton luminosities, etc.
