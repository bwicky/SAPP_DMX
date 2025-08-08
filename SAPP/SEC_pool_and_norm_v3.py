from opentrons import protocol_api
import numpy as np

metadata = {
    'protocolName': '### PROTOCOL_NAME ###',
    'author': 'Basile Wicky <bwicky@uw.edu>',
    'apiLevel':'2.5',
    'description': 'Pool SEC fractions and normalize concentrations (optional)',
}

#======================
# TRANSFER DEFINITIONS
#======================

instrument = '### INSTRUMENT ###' # {akta or hplc or cwby} instrument used for SEC -- determines the plate-type of the fractions. Use cwby if normalizing pre-pooled fractions.
normalize = ### NORMALIZE ### # whether or not to normalize concentrations of pooled fractions.
mixing = True # set to False to overwrite the automatic assignemnt of this variable.
debug = False # use for simulation purposes only.

"""
Format of transfer definitions (list of lists):
[
[['frac_1', 'frac_2', ...], [frac_1_vol, frac_2_vol, ...], 'destination_plate', 'destination_well', buffer_volume],
]

Format of frac_x strings:
'akta':'1.A.1'
'hplc':'1-P1-A1'
'cwby':'1_A1'
"""

transfers = np.array([
### TRANSFERS ###
], dtype=object)

def run(protocol: protocol_api.ProtocolContext):

    protocol.set_rail_lights(True)

    # Global defintions.
    z_off = 1.5 # mm -- tip offset from bottom of well.
    z_in = 3 # mm -- tip insertion below the liquid surface.

    # Labware definitions.
    s300 = protocol.load_instrument('p300_single_gen2', 'right')
    s300.flow_rate.aspirate = 300 # uL/s
    s300.flow_rate.dispense = 300 # uL/s
    s300.flow_rate.blow_out = 300 # uL/s
    
    if debug: # use debug plate definitions for simulation purposes only.
        plate_types = {
            'akta':'usascientific_96_wellplate_2.4ml_deep',
            'hplc':'corning_384_wellplate_112ul_flat',
            'cwby':'corning_96_wellplate_360ul_flat',
        }  
    
    else:
        plate_types = {
            'akta':'abgene_96_wellplate_2200ul',
            'hplc':'greiner_384_wellplate_250ul',
            'cwby':'nunc_96_wellplate_450ul',
        } 
    
    # Determine the number of destination plates.
    des_plates = np.unique(transfers[:, 2])
  
    # Determine the maximum volume of any well in the pooled plate
    if normalize:
        max_pool_vol = np.max([np.sum(t[1]) + t[4] for t in transfers])

    else: 
        max_pool_vol = np.max([np.sum(t[1]) for t in transfers])

    if max_pool_vol > 2000:
        protocol.pause('CAUTION!!! Some pooled wells will overflow!!!')

    # Determine the number of source plates.
    if instrument == 'akta':
        src_plates = list(sorted(np.unique([x.split('.')[0] for x in np.hstack(transfers[:, 0])])))
    
    elif instrument == 'hplc':
        src_plates = list(sorted(np.unique(['-'.join(x.split('-')[:2]) for x in np.hstack(transfers[:, 0])])))

    elif instrument == 'cwby':
        src_plates = list(sorted(np.unique([x.split('_')[0] for x in np.hstack(transfers[:, 0])])))

    else:
        print('Instrument not recognized.')
      
    # Assign plate-type and plate locations to deck.
    plate2deck = {0:10, 1:11, 2:7, 3:8, 4:4, 5:5}
    frac_ptype = plate_types[instrument]
    frac_plates = {}
    for i, p_name in enumerate(src_plates):  
        frac_plates[p_name] = protocol.load_labware(frac_ptype, plate2deck[i], label=f'Fractions {p_name}\n{frac_ptype}')
       
    tips = {}
    pool_ptype = plate_types['akta'] if max_pool_vol > 405 else plate_types['cwby']
    pool_plates = {}
    for i, p_name in enumerate(des_plates):
        tips[p_name] = protocol.load_labware('opentrons_96_tiprack_300ul', 9-(i*3), label=f'Tips {p_name}')
        pool_plates[p_name] = protocol.load_labware(pool_ptype, i+1, label=f'Destination {p_name}\n{pool_ptype}')
      

    def V2H(vol, z_insert, z_offset, plate):
        '''
        For 96-deepwell plate (AB0932); 39 x 8.4 x 8.4 mm
        For 384-deepwell plate (Greiner 781270); 19.3 x 3.8 x 3.8 mm
        '''
        h_factor = {
                'usascientific_96_wellplate_2.4ml_deep' :0.0177273,
                'corning_384_wellplate_112ul_flat'      :0.0804167,
                'corning_96_wellplate_360ul_flat'       :0.0,
                'abgene_96_wellplate_2200ul'            :0.0177273,
                'greiner_384_wellplate_250ul'           :0.0804167,
                'nunc_96_wellplate_450ul'               :0.0,

        }

        h = h_factor[plate] * vol # height of the liquid.
        h_corr = h - z_insert # height corrected for tip insertion.
        h_actual = h_corr if h_corr > z_offset else z_offset # ensure that min height is z_offset.
        return h_actual

    # Add buffer if normalizing concentrations.
    if normalize:
        buffer = protocol.load_labware('agilent_1_reservoir_290ml', 3, label='Buffer')
        s300.pick_up_tip(tips[des_plates[0]]['A1'])
        s300.aspirate(300, buffer['A1'].bottom(4))

        for t in transfers:
            _, _, des_p, des_w, buffer_vol = t

            if buffer_vol == 0:
                pass

            else:
            
                if buffer_vol <= s300.current_volume:
                    s300.dispense(buffer_vol, pool_plates[des_p][des_w].bottom(z_off))

                else:

                    if buffer_vol <= 300: 
                        s300.aspirate(300-s300.current_volume, buffer['A1'].bottom(4))
                        s300.dispense(buffer_vol, pool_plates[des_p][des_w].bottom(z_off))

                    else:
                        full_transfers = int(buffer_vol // 300)
                        remaining_vol = buffer_vol - (full_transfers * 300)
                        added_vol = 0
                        for _ in range(full_transfers):
                            s300.aspirate(300-s300.current_volume, buffer['A1'].bottom(4))
                            added_vol += 300
                            s300.dispense(300, pool_plates[des_p][des_w].bottom(V2H(added_vol, z_in, z_off, pool_ptype)))
                        
                        s300.aspirate(300, buffer['A1'].bottom(4))
                        added_vol += remaining_vol
                        s300.dispense(remaining_vol, pool_plates[des_p][des_w].bottom(V2H(added_vol, z_in, z_off, pool_ptype)))
        
        s300.blow_out(buffer['A1'])
        s300.return_tip()
    
    # Pool.   
    for t in transfers:
        fractions, frac_vol, des_p, des_w, buffer_vol = t
        tot_vol = np.sum(frac_vol) + buffer_vol
        s300.pick_up_tip(tips[des_p][des_w])
        added_vol = buffer_vol

        for frac, vol in list(zip(fractions, frac_vol)):

            if instrument == 'hplc':
                src_p, src_w = '-'.join(frac.split('-')[:2]), frac.split('-')[2]

            elif instrument == 'akta':
                src_p, src_w = frac.split('.')[0], ''.join(frac.split('.')[1:])

            elif instrument == 'cwby':
                src_p, src_w = frac.split('_')
            
            else:
                print('Instrument not recognized.')

            
            if instrument == 'hplc':
                # Step draw to avoid overflowing the well.
                # Only setup for taking everything from the well. 
                s300.aspirate(100, frac_plates[src_p][src_w].bottom(8.3)) 
                s300.aspirate(150, frac_plates[src_p][src_w].bottom(z_off))
                added_vol += vol
                s300.dispense(250, pool_plates[des_p][des_w].bottom(V2H(added_vol, z_in, z_off, pool_ptype)))
                
            else:
                full_transfers = int(vol // 300)
                remaining_vol = vol - (full_transfers * 300)
                for _ in range(full_transfers):
                    s300.aspirate(300, frac_plates[src_p][src_w].bottom(z_off))
                    added_vol += 300
                    s300.dispense(300, pool_plates[des_p][des_w].bottom(V2H(added_vol, z_in, z_off, pool_ptype)))
           
                s300.aspirate(remaining_vol, frac_plates[src_p][src_w].bottom(z_off))
                added_vol += remaining_vol
                s300.dispense(remaining_vol, pool_plates[des_p][des_w].bottom(V2H(added_vol, z_in, z_off, pool_ptype)))


                                                           
        if ((len(fractions) > 1) | ((buffer_vol != 0) & normalize)) & mixing: # mix if multiple fractions pooled or buffer was added.
            s300.mix(
                (lambda x: 3 if x <= 600 else (1.5 * tot_vol) // 300)(tot_vol), 
                (lambda x: x / 2 if x <= 600 else 300)(tot_vol), 
                pool_plates[des_p][des_w].bottom(V2H(tot_vol, z_in, z_off, pool_ptype))
            ) 

        s300.drop_tip()
    
    protocol.set_rail_lights(False)
from opentrons import protocol_api
