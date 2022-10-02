#!/usr/local/bin/vmd


proc rendering_high_res {file_name} {
	
	render Tachyon $file_name\.dat "~\Dropbox\Desktop\Desktop\Activation_potential\Potential_calc\force_field\VMD\Fig" -aasamples 10 -auto_skylight 0.7 -fullshade %s -format BMP -res 6400 3600 -o $file_name\.bmp
}

proc rendering_low_res {file_name} {

	render snapshot $file_name exlporer %s

}