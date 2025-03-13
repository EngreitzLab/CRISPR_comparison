MAX_MEM_MB = 250 * 1000  # 250GB

def determine_mem_mb(wildcards, input, attempt, min_gb=8):
	# Memory resource calculator for snakemake rules
	input_size_mb = input.size_mb
	if ".gz" in str(input):
		input_size_mb *= 8  # assume gz compressesed the file <= 8x
	attempt_multiplier = 2 ** (attempt - 1)  # Double memory for each retry
	mem_to_use_mb = attempt_multiplier *  max(4 * input_size_mb, min_gb * 1000)
	return min(mem_to_use_mb, MAX_MEM_MB)