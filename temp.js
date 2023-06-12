#!/usr/bin/env k8

var getopt = function(args, ostr) {
		var oli; // option letter list index
		if (typeof(getopt.place) == 'undefined')
				getopt.ind = 0, getopt.arg = null, getopt.place = -1;
		if (getopt.place == -1) { // update scanning pointer
					if (getopt.ind >= args.length || args[getopt.ind].charAt(getopt.place = 0) != '-') {
									getopt.place = -1;
									return null;
								}
					if (getopt.place + 1 < args[getopt.ind].length && args[getopt.ind].charAt(++getopt.place) == '-') { // found "--"
									++getopt.ind;
									getopt.place = -1;
									return null;
								}
				}
		var optopt = args[getopt.ind].charAt(getopt.place++); // character checked for validity
		if (optopt == ':' || (oli = ostr.indexOf(optopt)) < 0) {
					if (optopt == '-') return null; //  if the user didn't specify '-' as an option, assume it means null.
					if (getopt.place < 0) ++getopt.ind;
					return '?';
				}
		if (oli+1 >= ostr.length || ostr.charAt(++oli) != ':') { // don't need argument
					getopt.arg = null;
					if (getopt.place < 0 || getopt.place >= args[getopt.ind].length) ++getopt.ind, getopt.place = -1;
				} else { // need an argument
							if (getopt.place >= 0 && getopt.place < args[getopt.ind].length)
									getopt.arg = args[getopt.ind].substr(getopt.place);
							else if (args.length <= ++getopt.ind) { // no arg
											getopt.place = -1;
											if (ostr.length > 0 && ostr.charAt(0) == ':') return ':';
											return '?';
										} else getopt.arg = args[getopt.ind]; // white space
							getopt.place = -1;
							++getopt.ind;
						}
		return optopt;
}


Interval = {};

Interval.sort = function(a)
{
		if (typeof a[0] == 'number')
				a.sort(function(x, y) { return x - y });
		else a.sort(function(x, y) { return x[0] != y[0]? x[0] - y[0] : x[1] - y[1] });
}

Interval.merge = function(a, sorted)
{
		if (typeof sorted == 'undefined') sorted = true;
		if (!sorted) Interval.sort(a);
		var k = 0;
		for (var i = 1; i < a.length; ++i) {
					if (a[k][1] >= a[i][0])
							a[k][1] = a[k][1] > a[i][1]? a[k][1] : a[i][1];
					else a[++k] = a[i].slice(0);
				}
		a.length = k + 1;
}

Interval.index_end = function(a, sorted)
{
		if (a.length == 0) return;
		if (typeof sorted == 'undefined') sorted = true;
		if (!sorted) Interval.sort(a);
		a[0].push(0);
		var k = 0, k_en = a[0][1];
		for (var i = 1; i < a.length; ++i) {
					if (k_en <= a[i][0]) {
									for (++k; k < i; ++k)
											if (a[k][1] > a[i][0])
													break;
									k_en = a[k][1];
								}
					a[i].push(k);
				}
}

Interval.find_intv = function(a, x)
{
		var left = -1, right = a.length;
		if (typeof a[0] == 'number') {
					while (right - left > 1) {
									var mid = left + ((right - left) >> 1);
									if (a[mid] > x) right = mid;
									else if (a[mid] < x) left = mid;
									else return mid;
								}
				} else {
							while (right - left > 1) {
											var mid = left + ((right - left) >> 1);
											if (a[mid][0] > x) right = mid;
											else if (a[mid][0] < x) left = mid;
											else return mid;
										}
						}
		return left;
}


Interval.find_ovlp = function(a, st, en)
{
		if (a.length == 0 || st >= en) return [];
		var l = Interval.find_intv(a, st);
		var k = l < 0? 0 : a[l][a[l].length - 1];
		var b = [];
		for (var i = k; i < a.length; ++i) {
					if (a[i][0] >= en) break;
					else if (st < a[i][1])
							b.push(a[i]);
				}
		return b;
}


function paf_liftover(args)
{
		function read_bed(fn, to_merge)
		{
					if (fn == null) return null;
					if (typeof to_merge == 'undefined') to_merge = true;
					var file = fn == '-'? new File() : new File(fn);
					var buf = new Bytes();
					var bed = {};
					while (file.readline(buf) >= 0) {
									var t = buf.toString().split("\t");
									if (bed[t[0]] == null) bed[t[0]] = [];
									bed[t[0]].push([parseInt(t[1]), parseInt(t[2])]);
								}
					buf.destroy();
					file.close();

					for (var chr in bed) {
									Interval.sort(bed[chr]);
									if (to_merge)
											Interval.merge(bed[chr], true);
									Interval.index_end(bed[chr], true);
								}
					return bed;
				}

		var re_cigar = /(\d+)([MID])/g, re_tag = /^(\S\S):([AZif]):(\S+)$/;
		var c, to_merge = false, min_mapq = 5, min_len = 50000, max_div = 2.0;
		var re = /(\d+)([MID])/g;
		while ((c = getopt(args, "mq:l:d:")) != null) {
					if (c == 'm') to_merge = true;
					else if (c == 'q') min_mapq = parseInt(getopt.arg);
					else if (c == 'l') min_len = parseInt(getopt.arg);
					else if (c == 'd') max_div = parseFloat(getopt.arg);
				}
		if (args.length - getopt.ind < 2) {
					print("Usage: paftools.js liftover [options] <aln.paf> <query.bed>");
					print("Options:");
					print("  -q INT    min mapping quality [" + min_mapq + "]");
					print("  -l INT    min alignment length [" + min_len + "]");
					print("  -d FLOAT  max sequence divergence (>=1 to disable) [1]");
					exit(1);
				}
		var bed = read_bed(args[getopt.ind+1], to_merge);

		var file = args[getopt.ind] == '-'? new File() : new File(args[getopt.ind]);
		var buf = new Bytes();
		while (file.readline(buf) >= 0) {
					var t = buf.toString().split("\t");

					if (bed[t[0]] == null) continue; 

					
					var m, tp = null, cg = null;
					for (var i = 12; i < t.length; ++i) {
									if ((m = re_tag.exec(t[i])) != null) {
														if (m[1] == 'tp') tp = m[3];
														else if (m[1] == 'cg') cg = m[3];
													}
								}
					if (tp != 'P' && tp != 'I') continue; 
					if (cg == null) throw Error("unable to find the 'cg' tag");

					
					for (var i = 1; i <= 3; ++i)
							t[i] = parseInt(t[i]);
					for (var i = 6; i <= 11; ++i)
							t[i] = parseInt(t[i]);
					if (t[11] < min_mapq || t[10] < min_len) continue;
					var regs = Interval.find_ovlp(bed[t[0]], t[2], t[3]);
					if (regs.length == 0) continue; 
					if (max_div >= 0.0 && max_div < 1.0) {
									var n_gaps = 0, n_opens = 0;
									while ((m = re_cigar.exec(cg)) != null)
											if (m[2] == 'I' || m[2] == 'D')
													n_gaps += parseInt(m[1]), ++n_opens;
									var n_mm = t[10] - t[9] - n_gaps;
									var n_diff2 = n_mm + n_opens;
									if (n_diff2 / (n_diff2 + t[9]) > max_div)
											continue;
								}

					
					var a = [], r = [], strand = t[4];
					for (var i = 0; i < regs.length; ++i) {
									var s = regs[i][0], e = regs[i][1];
									if (strand == '+') {
														a.push([s, 0, i, -2]);
														a.push([e - 1, 1, i, -2]);
													} else {
																		a.push([t[1] - e, 0, i, -2]);
																		a.push([t[1] - s - 1, 1, i, -2]);
																	}
									r.push([-2, -2]);
								}
					a.sort(function(x, y) { return x[0] - y[0] });

					
					var k = 0, x = t[7], y = strand == '+'? t[2] : t[1] - t[3];
					while ((m = re_cigar.exec(cg)) != null) { 
									var len = parseInt(m[1]);
									if (m[2] == 'D') { 
														x += len;
														continue;
													}
									while (k < a.length && a[k][0] < y) ++k; 
									for (var i = k; i < a.length; ++i) {
														if (y <= a[i][0] && a[i][0] < y + len)
																a[i][3] = m[2] == 'M'? x + (a[i][0] - y) : x;
														else break;
													}
									y += len;
									if (m[2] == 'M') x += len;
								}
					if (x != t[8] || (strand == '+' && y != t[3]) || (strand == '-' && y != t[1] - t[2]))
							throw Error("CIGAR is inconsistent with mapping coordinates");

					
					for (var i = 0; i < a.length; ++i) {
									if (a[i][1] == 0) r[a[i][2]][0] = a[i][3];
									else r[a[i][2]][1] = a[i][3] + 1; 
								}
					for (var i = 0; i < r.length; ++i) {
									var name = [t[0], regs[i][0], regs[i][1]].join("_");
									if (r[i][0] < 0) name += "_t5", r[i][0] = t[7];
									if (r[i][1] < 0) name += "_t3", r[i][1] = t[8];
									print(t[5], r[i][0], r[i][1], name, 0, strand,cg, t[7]);
								}
				}
		buf.destroy();
		file.close();
}

function main(args)
{
		var cmd = args.shift();
		if (cmd == 'liftover' || cmd == 'liftOver') paf_liftover(args);
		else throw Error("unrecognized command: " + cmd);
}

main(arguments);
