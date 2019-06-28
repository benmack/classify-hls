from nasa_hls import hls_qa_layer_to_mask


path__qa = snakemake.input[0]
path__clear = snakemake.output[0]

# valid ids
valid = [  0,   4,  16,  20,  32,  36,  48,  52,  64,  68,  80,  84,  96,
         100, 112, 116, 128, 132, 144, 148, 160, 164, 176, 180, 192, 196,
         208, 212, 224, 228, 240, 244]

clear = hls_qa_layer_to_mask(qa_path=path__qa,
                             qa_valid=valid,
                             keep_255=True,
                             mask_path=path__clear,
                             overwrite=False)

assert clear == 0