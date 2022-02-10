def get_site_positions(predicted_target_info_data):
    UTR_positions = []
    for index, row in predicted_target_info_data.iterrows():
        UTR_start = row["UTR start"]
        UTR_end = row["UTR end"]
        UTR_start_end = f"{UTR_start}:{UTR_end}"
        UTR_positions.append(UTR_start_end)

    return UTR_positions
