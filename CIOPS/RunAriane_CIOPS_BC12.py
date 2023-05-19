"""
Loop over dates, setting up a series of forcing links, running Ariane starting particles over multiple days
and saving the statistic results"""

import argparse
import arrow
import os
import subprocess
import datetime as dt
import pytz

TARGET_TMPL = 'CIOPS_1h_{:06d}_grid_{:s}.nc'
FILENAME_TMPL = 'BC12_1h_grid_{:s}_2D_{:%Y%m%d}_{:%Y%m%d}.nc'
FILENAME_TMPL2 = '{:s}_{:%Y%m%d}.nc'
SUBDIR_TMPL = '{:%Y%m%d}00'
SUBDIR_TMPL2 = '{:%Y%m}'
NEW_SUBDIR_TMPL = 'back_sahubdy_{:%d%b%y}'

utc=pytz.UTC
folders = [(dt.datetime(2015,11,22)+dt.timedelta(days=7*(i+1))).replace(tzinfo=utc) for i in range(int(214))]

def make_links(rundate, runlength):
    tardir = 'Links'
    for grid in ['T', 'U_new', 'V_new', 'S_new', 'T_new']:
        for fileno in range(runlength):
            target = TARGET_TMPL.format(fileno+1, grid)
            date = rundate.shift(days=+fileno).datetime
            if grid == 'T':
                dir = '/ocean/mdunphy/CIOPSW-BC12/'
                for i in range(len(folders)-1):
                    d = date.replace(tzinfo=utc)
                    if d >= folders[i] and d < folders[i+1]:
                        folderdate = folders[i+1]
                    elif d < folders[0]:
                        folderdate = folders[0]
                link = FILENAME_TMPL.format(grid, date, date)
                subdir = SUBDIR_TMPL.format(folderdate)
            else:
                dir = '/ocean/rbeutel/data/'
                link = FILENAME_TMPL2.format(grid, date)
                subdir = SUBDIR_TMPL2.format(date)
            try:
                os.unlink(os.path.join(tardir, target))
            except FileNotFoundError:
                pass
            os.symlink(os.path.join(dir, subdir, link),
                       os.path.join(tardir, target))


def run_ariane():
    with open('babypoo', 'wt') as stdout:
        with open('errpoo', 'wt') as stderr:
            subprocess.run(
                "/ocean/rbeutel/MOAD/ariane-2.3.0_03/bin/ariane",
                stdout=stdout, stderr=stderr,
                universal_newlines=True)


def rename_results(rundate=arrow.utcnow(), subdir='',
                   nday=1, labeltype='date'):
    finaldir = '/ocean/rbeutel/MOAD/analysis-becca/Ariane/'
    filelist = ['ariane_memory.log', 'ariane_statistics_quantitative.nc',
                'final_pos.txt', 'init_pos.txt', 'output',
                'ariane_positions_quantitative.nc', 'final.sav',
                'init.sav', 'stats.txt']
    for filename in filelist:
        if labeltype == 'date':
            newdir = NEW_SUBDIR_TMPL.format(rundate.datetime).lower()
        else:
            newdir = ''
        print(os.path.join(finaldir, subdir, newdir, filename))
        if not os.path.exists(os.path.join(finaldir, subdir, newdir)):
                os.makedirs(os.path.join(finaldir, subdir, newdir))
        os.rename(filename, os.path.join(finaldir, subdir, newdir, filename))


def main(args):
    initialrundate = arrow.get(args.initialrundate, 'YYYY-MM-DD')
    for nday in range(args.numberofdays):
        rundate = initialrundate.shift(days=+nday)
        print ('Start', rundate)
        if args.forback == 'forward':
            startfile = rundate
        elif args.forback == 'backward':
            startfile = rundate.shift(days=-(args.runlength+args.releaseparticledays-1))
        make_links(startfile, args.runlength+args.releaseparticledays)
        print ('Startfile', startfile)
        run_ariane()
        rename_results(rundate=rundate, subdir=args.runtype)
        print ('End', rundate)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'initialrundate', help='Date to start from as YYYY-MM-DD', type=str)
    parser.add_argument(
        'numberofdays', help='Number of different dates to do', type=int)
    parser.add_argument(
        'runlength', help='Number of days to track particles', type=int)
    parser.add_argument(
        'releaseparticledays', help='Number of days over which to release particles for one run', type=int)
    parser.add_argument('forback', help='Run forward or backward', type=str)
    parser.add_argument('runtype',
                        help='FluxesSouth FluxesNorth BackSouth BackNorth',
                        type=str)
    args = parser.parse_args()
    main(args)
