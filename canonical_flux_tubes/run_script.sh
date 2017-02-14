#! /bin/bash

# Temperature Density Tubes
timeout 17m smt run --reason "Look at temperature and density from default view." --density_tubes --temperature_tubes --view default --output_prefix temperature_density_tubes_ 2016-10-25-20-10
timeout 17m smt run ---reason "Look at temperature and density from default lower angle view." --density_tubes --temperature_tubes --view default_lower_angle --output_prefix temperature_density_tubes_ 2016-10-25-20-10
timeout 17m smt run --reason "Look at temperature and density along positive z view." --density_tubes --temperature_tubes --view positive_z --output_prefix temperature_density_tubes_ 2016-10-25-20-10
timeout 17m smt run --reason "Look at temperature and density along negative z view." --density_tubes --temperature_tubes --view negative_z --output_prefix temperature_density_tubes_ 2016-10-25-20-10
timeout 17m smt run --reason "Look at temperature and density along positive x view." --density_tubes --temperature_tubes --view positive_x --output_prefix temperature_density_tubes_ 2016-10-25-20-10

# Magnetic Density Tubes
timeout 17m smt run --reason "Look at magnetic and density from default view." --density_tubes --electron --view default --output_prefix electron_and_density_tubes_ 2016-10-25-20-10
timeout 17m smt run ---reason "Look at magnetic and density from default lower angle view." --density_tubes --electron --view default_lower_angle --output_prefix electron_and_density_tubes_ 2016-10-25-20-10
timeout 17m smt run --reason "Look at magnetic and density along positive z view." --density_tubes --electron --view positive_z --output_prefix electron_and_density_tubes_ 2016-10-25-20-10
timeout 17m smt run --reason "Look at magnetic and density along negative z view." --density_tubes --electron --view negative_z --output_prefix electron_and_density_tubes_ 2016-10-25-20-10
timeout 17m smt run --reason "Look at magnetic and density along positive x view." --density_tubes --electron --view positive_x --output_prefix electron_and_density_tubes_ 2016-10-25-20-10

# Magnetic Temperature Tubes
timeout 17m smt run --reason "Look at magnetic and temperature from default view." --temperature_tubes --electron --view default --output_prefix temperature_and_electron_tubes_ 2016-10-25-20-10
timeout 17m smt run ---reason "Look at magnetic and temperature from default lower angle view." --temperature_tubes --electron --view default_lower_angle --output_prefix temperature_and_electron_tubes_ 2016-10-25-20-10
timeout 17m smt run --reason "Look at magnetic and temperature along positive z view." --temperature_tubes --electron --view positive_z --output_prefix temperature_and_electron_tubes_ 2016-10-25-20-10
timeout 17m smt run --reason "Look at magnetic and temperature along negative z view." --temperature_tubes --electron --view negative_z --output_prefix temperature_and_electron_tubes_ 2016-10-25-20-10
timeout 17m smt run --reason "Look at magnetic and temperature along positive x view." --temperature_tubes --electron --view positive_x --output_prefix temperature_and_electron_tubes_ 2016-10-25-20-10

# ion Omega_i_raw_plus with for loop
for time_point in {0..249}
    do timeout 17m smt run --reason "Look at ion tubes Omega_i_raw_plus with varying start time." --ion --view default --start_time_point $time_point --output_prefix ion_tubes_ 2016-10-25-20-10
done

# ion Omega_i_raw_plus_density_dependence for loop
for time_point in {0..249}
    do timeout 17m smt run --reason "Look at ion tubes Omega_i_raw_plus_density_dependence with varying start time." --ion --view default --start_time_point $time_point --omega_to_use Omega_i_raw_plus_density_dependence --output_prefix ion_tubes_ 2016-10-25-20-10
done

# ion Omega_i_raw_plus_times_density for loop
for time_point in {0..249}
    do timeout 17m smt run --reason "Look at ion tubes Omega_i_raw_plus_density with varying start time." --ion --view default --start_time_point $time_point --omega_to_use Omega_i_raw_plus_density --output_prefix ion_tubes_ 2016-10-25-20-10
done
