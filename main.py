import argparse
import json


def main():

    # Parse information from the command line
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("namelist")
    args = parser.parse_args()

    file_namelist = open(args.namelist).read()
    namelist = json.loads(file_namelist)
    del file_namelist

    main_run(namelist)

    return


def main_run(case_ic):
    import Simulation

    simulation = Simulation.Simulation()
    simulation.initialize(case_ic)
    simulation.run()

    return


if __name__ == "__main__":
    main()
