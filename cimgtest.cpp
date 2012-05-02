/*
 * cimgtest.cpp
 *
 *  Created on: Apr 27, 2012
 *      Author: tobias
 */

#include <algorithm>
#include <vector>
#include <map>
#include <string>
#include "CImg.h"

using namespace cimg_library;

float fabsf(float);

void cluster_search(const CImg<float> & fimage, CImg<unsigned char> & cluster,
		int x, int y, float sameness_limit)
{
	std::vector<std::pair<int,int> > open_set(1, std::pair<int,int>(x,y));
	while (open_set.size()) {
		x = open_set.back().first;
		y = open_set.back().second;
		open_set.pop_back();
		cluster(x,y) = 255;
		for (int xx = x-1; xx <= x+1; ++xx) for (int yy = y-1; yy <= y+1; ++yy) {
			if (xx >= 0 && xx < cluster.width() && yy >= 0 && yy < cluster.height())
				if (cluster(xx,yy) == 0 && fabsf(fimage(x,y) - fimage(xx,yy)) < sameness_limit) {
					cluster(xx, yy) = 128;
					open_set.push_back(std::pair<int,int>(xx,yy));
		        }
		}
	}
}

CImg<float> rgb_to_gray(const CImg<unsigned char> & img) {
	int w = img.width();
	int h = img.height();
	CImg<float> fimage(w,h,1,1,0.0f);
	for (int x = 0; x < w; ++x) for (int y = 0; y < h; ++y)
		fimage(x,y) =
				(img(x,y,0,0) * 0.30f + img(x,y,0,1) * 0.59f + img(x,y,0,2) * 0.11f)
				/ 256.0f;
	return fimage;
}



int max_hist_bin_index(const CImg<float> & hist)
{
    float maxval = -1e6;
    int maxind = -1;
    for(int x = 0;x < hist.width();++x){
        if(hist(x, 0) > maxval){
            maxval = hist(x, 0);
            maxind = x;
        }
    }
    return maxind;
}

int center_square_max_hist_bin(const CImg<float> & fimage, const int bins) {
	int w = fimage.width();
	int h = fimage.height();

	int lsquare = std::min(w,h) / 3;
	int left = (w - lsquare) / 2;
	int top = (h - lsquare) / 2;
	CImg<float> square = fimage.get_crop(left,top,left+lsquare,top+lsquare);

	CImg<float> hist = square.get_histogram(bins, 0.0f, 1.0f);
    return max_hist_bin_index(hist);
}

CImg<unsigned char> find_biggest_cluster(int w, int h, const float initial_white, CImg<float> fimage, const float initial_tolerance, const float sameness_limit)
{
    CImg<unsigned char> clustered_flag(w, h, 1, 1, 0), biggest_cluster(clustered_flag), current_cluster(clustered_flag);
    double biggest_cluster_size = 0;
    for(int x = 0;x < w;++x)
        for(int y = 0;y < h;++y){
            if(clustered_flag(x, y) == 0 && fabsf(initial_white - fimage(x, y)) < initial_tolerance){
                current_cluster.fill((unsigned char)(0));
                cluster_search(fimage, current_cluster, x, y, sameness_limit);
                double current_cluster_size = current_cluster.sum() / 255;
                clustered_flag |= current_cluster;
                if(current_cluster_size > biggest_cluster_size){
                    printf("replaces biggest cluster\n");
                    biggest_cluster.assign(current_cluster);
                    biggest_cluster_size = current_cluster_size;
                }
            }

        }
    return biggest_cluster;
}

void find_white_area(const CImg<float> & fimage)
{
    int w = fimage.width();
    int h = fimage.height();
    const int bins = 64;
    const float initial_tolerance = 1.0f / bins / 2;
    const float initial_white =
    		1.0f / bins * center_square_max_hist_bin(fimage, bins) + initial_tolerance;
    const float sameness_limit = 1.0f / 256 * 8; //4;
    CImg<unsigned char> biggest_cluster =
    		find_biggest_cluster(w, h, initial_white, fimage, initial_tolerance, sameness_limit);
    biggest_cluster.save_png("/tmp/x.png", 1);
    biggest_cluster.display("White Area");
}

int max_hist_bin_width(const CImg<float> & hist, int maxind)
{
    float maxval = hist(maxind, 0);
    int width = 2, searchdownind = maxind - 1, searchupind = maxind + 1;
    while(searchdownind > 0 && hist(searchdownind, 0) > maxval / 2){
        --searchdownind;
        ++width;
    }
    while(searchupind < hist.width() && hist(searchupind, 0) > maxval / 2){
        ++searchupind;
        ++width;
    }
    return width;
}

void crop_at_white_margin(CImg<float> & fimage, const float white = 0.67)
{
	int w = fimage.width(), h = fimage.height();
	enum {north,south,east,west};
	int crop[4];
	for (int blowing = north; blowing <= west; ++blowing) {
		printf("searching %d\n", blowing);
		crop[blowing] = 0;
		int max_search, max_dark_edges;
		bool found_white = false;
		if (blowing <= south) {
			max_search = h / 5;
			max_dark_edges = w / 15;

		} else {
			max_search = w / 5;
			max_dark_edges = h / 15;
		}
		for (int search = 0; search < max_search; ++search) {
			CImg<float> test;
			switch(blowing) {
			case north:
				test = fimage.get_crop(0,search,w-1,search);
				break;
			case south:
				test = fimage.get_crop(0,h-search-1,w-1,h-search-1);
				break;
			case west:
				test = fimage.get_crop(search,0,search,h-1).transpose();
				break;
			case east:
				test = fimage.get_crop(w-search-1,0,w-search-1,h-1).transpose();
				break;
			}
			// count dark pixels, ignore edge case as defined
			int low_edge_last_dark = -1, high_edge_last_dark = -1;
			for (int edge_ind = 0; edge_ind < max_dark_edges; ++edge_ind) {
				if (test(edge_ind) < white)
					low_edge_last_dark = edge_ind;
				if (test(test.width()-edge_ind-1) < white)
					high_edge_last_dark = edge_ind;
			}
			bool dark = low_edge_last_dark + high_edge_last_dark + 2 > max_dark_edges;
			for (int ind = max_dark_edges; ind < test.width() - max_dark_edges; ++ind)
				dark = dark || (test(ind) < white);
			if (found_white) {
				crop[blowing] = search;
				if (dark) break;
			}
			found_white = !dark;
		}
	}
	fimage.crop(crop[west],crop[north],w-crop[east]-1,h-crop[south]-1);
}

void white_statistics(CImg<float> & fimage, const int bins,
		const int max_adjacent_white_diff = 2,
		const int M = 11, const int N = 15)
{
    int w = fimage.width();
    int h = fimage.height();

	CImg<float> hist = fimage.get_histogram(bins, 0.0f, 1.0f);
	int gmaxind = max_hist_bin_index(hist);
    int gwidth = max_hist_bin_width(hist, gmaxind);
    int min_white = std::max(1, gmaxind - gwidth);
    int max_white = std::min(bins - 1, gmaxind + gwidth);

	printf("white= %d +- %d\n", gmaxind, gwidth);

	int tile_w = w / N;
	int tile_h = h / M;
	CImg<int> tiles_white(N,M,1,1,-1);
	CImg<int> tiles_white_width(tiles_white);

	for (int m = 0; m < M; ++m) {
		for (int n = 0; n < N; ++n) {
			CImg<float> hist = fimage.get_crop(n * tile_w, m * tile_h,
					(n+1)*tile_w-1, (m+1)*tile_h-1).histogram(bins,0.0f,1.0f);
			tiles_white(n,m) = max_hist_bin_index(hist);
			tiles_white_width(n,m) = max_hist_bin_width(hist, tiles_white(n,m));

			printf("%d ",tiles_white(n,m));
		}
		printf("\n");
	}

	// white value correction;
	// Determine values that need not be corrected.
	std::vector<std::pair<int,int> > open_set, correct_these_set;
	CImg<unsigned char> tiles_checked(N,M,1,1,0);
	for (int m = 0; m < M; ++m)	for (int n = 0; n < N; ++n) {
		if (tiles_white(n,m) == gmaxind) {
			open_set.push_back(std::pair<int,int>(n,m));
			tiles_checked(n,m) = 128;
		}
	}
	if (open_set.size() == 0) {
		fprintf(stderr, "Empty open set\n");
		exit(1);
	}
	while (open_set.size()) {
		int n = open_set.front().first;
		int m = open_set.front().second;
		open_set.erase(open_set.begin());
		for (int nn = n-1; nn <= n+1; ++nn) for (int mm = m-1; mm<=m+1; ++mm) {
			if (nn >= 0 && nn < N && mm >= 0 && mm < M)
				if (tiles_checked(nn,mm) == 0) {
					if (tiles_white(nn,mm) >= min_white && tiles_white(nn,mm) <= max_white &&
						abs(tiles_white(nn,mm) - tiles_white(n,m)) <= max_adjacent_white_diff) {
							open_set.push_back(std::pair<int,int>(nn,mm));
							tiles_checked(nn,mm) = 128;
					} else {
						correct_these_set.push_back(std::pair<int,int>(nn,mm));
					}
				}
		}
	}
	// Correct the rest
	for (unsigned k = 0; k < correct_these_set.size(); ++k) {
		int n = correct_these_set[k].first;
		int m = correct_these_set[k].second;

		if (tiles_checked(n,m) == 0) {
			int numerator = 0, denominator = 0;
			for (int nn = n-1; nn <= n+1; ++nn) for (int mm = m-1; mm<=m+1; ++mm)
				if (nn >= 0 && nn < N && mm >= 0 && mm < M) {
					if (tiles_checked(nn,mm)) {
						numerator += tiles_white(nn,mm);
						++denominator;
					}
					else
						correct_these_set.push_back(std::pair<int,int>(nn,mm));
				}
			tiles_white(n,m) = denominator ? (numerator / denominator) : gmaxind;
			tiles_checked(n,m) = 64;
		}
	}
	printf("corrected:\n");
	for (int m = 0; m < M; ++m) {
		for (int n = 0; n < N; ++n) {
			printf("%d ",tiles_white(n,m));
		}
		printf("\n");
	}

	// whiten
	CImg<bool> outside_white_corridor(w,h,1,1,false);
	for (int x = 0; x < w; ++x) for (int y = 0; y < h; ++y) {
		int m = std::min(y / tile_h, M-1);
		int n = std::min(x / tile_w, N-1);
		int white = tiles_white(n,m);
		int width = std::min(tiles_white_width(n,m), 3);

		if (fabs(white - fimage(x,y) * bins) > width)
			outside_white_corridor(x,y) = true;
		float white_lower_limit = float(white - width) / bins;
		fimage(x,y) = std::min(255.0/256,pow(fimage(x,y) / white_lower_limit, 6));
	}

	// crop 1: discard columns/lines with mostly non-white pixels
	int top,bot,lef,rig;
	for (top = 0;
			top < h && outside_white_corridor.get_crop(0,top,w-1,top).sum() > (w / 2.0);
			++top);
	for (bot = h-1;
			bot > top &&  outside_white_corridor.get_crop(0,bot,w-1,bot).sum() > (w / 2.0);
			--bot);
	for (lef = 0;
			lef < w && outside_white_corridor.get_crop(lef,0,lef,h-1).sum() > (h / 2.0);
				++lef);
	for (rig = w-1;
				rig > lef &&  outside_white_corridor.get_crop(rig,0,rig,h-1).sum() > (h / 2.0);
				--rig);
	fimage.crop(lef,top,rig,bot);

	// crop 2: search for white-ish margin
	crop_at_white_margin(fimage);

}
void black_statistics(const CImg<float> & fimage, const int bins,
		const int M = 11, const int N = 15)
{
    int w = fimage.width();
    int h = fimage.height();

	CImg<float> hist = fimage.get_histogram(bins, 0.0f, 1.0f);
	int gmaxind = max_hist_bin_index(hist);
    int gwidth = max_hist_bin_width(hist, gmaxind);

	printf("white= %d +- %d\n", gmaxind, gwidth);

	int tile_w = w / N;
	int tile_h = h / M;
	CImg<int> tiles(N,M,1,1,-1);
	for (int m = 0; m < M; ++m) {
		for (int n = 0; n < N; ++n) {
			CImg<float> hist = fimage.get_crop(n * tile_w, m * tile_h,
					(n+1)*tile_w-1, (m+1)*tile_h-1).histogram(bins,0.0f,1.0f);
			int maxind = max_hist_bin_index(hist);
			int width = max_hist_bin_width(hist, maxind);
			printf("%d+-%d ",maxind,width);
			if (maxind >= (gmaxind - gwidth) && maxind <= (gmaxind + gwidth))
				tiles(n,m) = maxind;

		}
		printf("\n");
	}
	// shepard interpolation;
	CImg<float> ftiles(tiles);
	for (int m = 0; m < M; ++m)
		for (int n = 0; n < N; ++n) {
			if (tiles(n,m) < 0) {
				float weightsum = 0;
				ftiles(n,m) = 0;
				for (int x = 0; x < N; ++x) for (int y = 0; y < M; ++y)
					if (tiles(x,y) >= 0) {
						float weight = 1 / sqrt((x-n)*(x-n) + (y-m)*(y-m));
						weightsum += weight;
						ftiles(n,m) += weight * tiles(x,y);
					}
				ftiles(n,m) /= weightsum;
			}
			else ftiles(n,m) = tiles(n,m);
		}

}


int main(int argc, char ** argv)
{
	std::string filename = "/home/tobias/Buch/Softwareentwicklung/Android/Android_Web_Apps/DSCN0477.JPG";
	if (argc >= 3)
		filename = argv[2];

	std::string mode = "find_white";
	if (argc >= 2)
		mode = argv[1];
	if (mode == "fine_white") {
		CImg<unsigned char> image(filename.c_str());
	    CImg<float> fimage = rgb_to_gray(image);
		find_white_area(fimage);
	} else if (mode == "white_statistics") {
		CImg<unsigned char> image(filename.c_str());
	    CImg<float> fimage = rgb_to_gray(image);

		const int bins = 64;
		white_statistics(fimage, bins);
		for (int fnameind = 3; fnameind < argc; ++fnameind) {
			CImg<unsigned char> image(argv[fnameind]);
			CImg<float> fimage = rgb_to_gray(image);
			white_statistics(fimage, bins);
		}
	}
	else if (mode == "black_statistics") {
		CImg<unsigned char> image(filename.c_str());
	    CImg<float> fimage = rgb_to_gray(image);

	    const int bins = 64;
		white_statistics(fimage, bins);
		for (int fnameind = 3; fnameind < argc; ++fnameind) {
			CImg<unsigned char> image(argv[fnameind]);
			CImg<float> fimage = rgb_to_gray(image);
			white_statistics(fimage, bins);
		}
	} else if (mode == "book") {
		std::string output_dir = argv[2];
		std::map<double, std::string> page_map;
		for (int i = 3; i < (argc - 1); i+= 2) {
			filename = argv[i];
			double page = atof(argv[i+1]);
			page_map[page] = filename;
		}
		for (std::map<double,std::string>::const_iterator i = page_map.begin();
				i != page_map.end(); ++i) {
			CImg<unsigned char> image((*i).second.c_str());
			CImg<float> fimage = rgb_to_gray(image);
			const int bins = 64;
			white_statistics(fimage, bins);
			char output_name[100] = "";
			double page = (*i).first;
			snprintf(output_name, 99, "/%07.2f.png", page);
			fimage *= 256;
			CImg<float> rotated(fimage.height(), fimage.width(), 1, 1);
			bool odd = false;
			if (page >= 1.0) odd = int(page) & 1;
			else odd = int(page * 100) & 1;
			for (int x = 0; x < fimage.width(); ++x) for (int y = 0; y < fimage.height(); ++y) {
				if (odd)
					rotated(rotated.width() - y - 1,x) = fimage(x,y);
				else
					rotated(y, rotated.height() - x - 1) = fimage(x,y);
			}
			rotated.save_png((output_dir + output_name).c_str(), 1);
		}
	}
	return 0;
}
