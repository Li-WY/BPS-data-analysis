.. _bps-protocol:

BPS Protocol
################################


Here's a protocol for doing Bacterial Positioning System, here about the
wetbench work. Consult the paper for more complete context.

.. contents::
    :depth: 2
    :local:
    :backlinks: top

Designing the plasmids and strains
==================================================================

Presently, our system requires elements to be present.
We are working on making it more general, but that is work in progress.

(diagrams and such)


This protocol assumes that cloning has been performed in recipient plasmids and
recipient cells that are compatible with the arrayed mating platform.

Recipient plasmids must contain a swapping cassette with I-SceI digestion 
sites and homology regions surrounding a counter-selectable cassette.

We have found that relE works effectively as a counter-selectable cassette.

We typically include a selectable marker (Hyg/clonNat) on the swapping 
cassette in the donor plasmid
to improve selection and for downstream applications, but this is not required.
Recipient cells must include an inducible lambda red and an inducible I-SceI
restriction enzyme on either a helper plasmid or on the cell genome.



Equipment and/or Consumables
====================================================================

* colony replicator 

    This key tool can be as simple as a toothpick, as complicated as a pinning 
    robot (such as a Singer ROTOR), or some where in between - like a 
    metal 96-well pinner, a carefully used multichannel pipette, or even a
    skillfully-used velvet.

    BPS requires physically combining colonies in known ways, so whatever tool
    accomplishes this for you best is what will largely determine the scale.

* Oxford Nanopore MinION sequencer, flowcell, and library preparation reagents

    BPS works best with long-read sequencing, so an ONT sequencer is nice to
    use. However, nowdays researchers can outsource their sequencing to 
    small companies such as Plasmidsaurus or Primordium.

    For reference, if you are considering doing your own sequencing then
    a MinION Flow Cell (R10.4.1) should cost about $900 each and be sufficient
    for ~10,000 plasmids. Library preparation should cost about $100,
    with additional costs for sample multiplexing of about $3.50/sample for 
    96 plex to $14/sample for 24plex. Thus to make maximal use of the ~750
    recipient barcoder strains here, this would require about 12 libraries
    multiplexed into one library prep - thus this depends on if you buy the
    24 or 96 -plex kit.

* colony picking

    If you are arraying colonies from a mixture, you will of course have to 
    pick these colonies into an array. To do this, you can use something as
    simple as toothpicks or pipette tips or as complicated as a dedicated
    automated colony picking robot. You may also desire an automated solution
    for rearraying or cherry-picking desired colonies from your BPS'd
    array for later use.

    For more see :ref:`tips about increasing scale <scale-tips>`.

* array agar plates

    Your arrays and the BPS strains can be isolated and handled on 
    agar pads or in multiwell plates.
    Agar pads/plates can be poured in Singer plus plates, or Omnitrays. 
    This method is tested with the following types of media, made using
    conventional technique (mixing dry ingredients in appropriate volume
    of deionized water, including 20g/L agar for solid media, 
    and autoclaving for ~20min before pouring in trays 
    or storing for later dispensing).

    * Donor plates: LB + Kanamycin + Hyg/clonNat
    * Recipient plates: LB + Sp + Gm + 2% Glucose
    * Mating plates: LB + Ara + IPTG
    * Selection plates: LB + Ara + Rha + Gm + Hyg/clonNat

General protocol
=====================================================

.. 
    This protocol should take X-Y days for Z plates of clones, at a total cost of 
    $A.BC for consumables.

Prepare strains for mating. 
-----------------------------------------------------

1. Prepare strains for mating. 
    Isolate/obtain the strains that you wish to "position" with BPS.

Colony picking and pre-growth (~20 minutes hands-on-time per 96 colonies, 6-24
hours cell growth) 

Thaw 96-well plate(s) containing glycerol stocks of donor
barcode strains by placing them on the bench at RT.


Let plates thaw while picking colonies.


If performing the procedure frequently, you can keep a stock
at 4ÂºC in liquid in 96-well plates or on an agar pad.


Pick single colonies
into 100 Âµl of LB with appropriate antibiotics in a 96-well plate.


Tip from
Octant: Pick colonies using 200 ÂµL pipette tips, put them back into the box,
and then use a multichannel pipette to inoculate the plate.


Replica plate the
thawed barcode plate(s) into 96-well plates using a 96-pin manual replicator.
Incubate donor and recipient cell arrays at 37ÂºC at 250 rpm for a minimum of 6
hours.


Warm mating plates for the next step at 37ÂºC.


[generate glycerol stock
for clones with DNA of interest] 


Mate strains and induce recombination
-----------------------------------------------------

Mating (~5 minutes hands-on-time, 4-6 hours
cell growth) Pin donor and recipient cells onto the same mating plate using a
96-pin manual replicator.


Replicating tools with guided pinning are recommended
to assure that donor cells and recipient cells make contact at each location.
Pin at densities of up to 384 position/plate for replicators that enable
pinning in this format.


For the quickest protocol, mating can be performed the
same day as pre-growth, but this can be a long day.


Alternatively, mating can
be performed the next morning.


IMPORTANT: Each location on a plate must have a
unique barcode â if pinning onto 384-position plates, you must have a minimum
of 384 barcodes.


 Incubate a minimum of 4 hours at 37ÂºC.


Warm selection plates
for the next step at 37ÂºC.


Selection of recombinants
-----------------------------------------------------

Selection (5 minutes hands-on-time, 12-24 hours

cell growth) Pin mating cells onto selection plates using a replicator.
Incubate overnight at 37ÂºC for 12-24 hours.



Collect samples
-----------------------------------------------------

Library prep (1.5-2 hours hands-on-time)
Collect colonies from plates by adding 5ml of liquid LB and scraping cells with a cell spreader off the agar surface
into the liquid LB. 

Resuspend cells into the LB using an automated pipettor and
add to the appropriate sized conical tube. 

Selection plates that contain
different barcode sets can be pooled together. 

Larger pools are more cost and
time efficient for subsequent steps.  

Vortex to mix pools well.  

Perform a
miniprep on each pool according to the kit protocol. 

Only 1.5-3 ml of pooled cells are required for each miniprep.  


Measure DNA concentration using your
preferred method.  

Subject plasmid pools to long-read sequencing
-----------------------------------------------------

outsource or diy

Prepare Oxford Nanopore libraries using the ONT rapid
barcoding kit according to the kit protocol. 

If you want to multiplex several
minipreps, use a unique ONT barcode for each miniprep and be sure to note which
ONT barcode marks each plate or set of plates. 

ONT barcodes (plate barcodes)
will later be used to demultiplex plates. 

Positional barcodes introduced by
arrayed conjugation can be used repeatedly so long as they are not repeated
within the same ONP barcode pool.

Sequencing

.. _scale-tips:

Tips and tricks for increasing scale of BPS
====================================================================


