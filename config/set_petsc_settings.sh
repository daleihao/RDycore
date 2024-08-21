#!/bin/sh

mach=
config=1

display_help() {
    echo "Usage: $0 " >&2
    echo
    echo "   -h, --help                        Display this message"
    echo "   --mach <pm-cpu|pm-gpu|frontier>   Supported machine name"
    echo "   --config <1|2|3>                  Configuration (1 = default)"
    echo
    echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    echo "Supported PETSc configurations"
    echo
    echo "For Perlmutter CPU (--mach pm-cpu): "
    echo "  --config 1: Without debugging, 64bit indices, and HDF5 1.14.3"
    echo "  --config 2: With debugging, 64bit indices, and HDF5 1.14.3"
    echo "  --config 3: Without debugging, 32bit indices, and HDF5 1.12.2.3"
    echo
    echo "For Perlmutter CPU (--mach pm-gpu): "
    echo "  --config 1: Without debugging, 64bit indices, and HDF5 1.14.3"
    echo "  --config 2: With debugging, 64bit indices, and HDF5 1.14.3"
    echo "  --config 3: Without debugging, 32bit indices, and HDF5 1.12.2.3"
    echo
    echo "For Frontier (--mach frontier): "
    echo "  --config 1: Without debugging, 64bit indices, and HDF5 1.14.3"
    echo "  --config 2: With debugging, 64bit indices, and HDF5 1.14.3"
    echo "  --config 3: Without debugging, 32bit indices, and HDF5 1.12.2.1"
    echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    echo

    return 1
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

while [ $# -gt 0 ]
do
  case "$1" in
    --mach ) mach="$2"; shift ;;
    --config) config="$2"; shift ;;
    -*)
      display_help
      exit 0
      ;;
    -h | --help)
      display_help
      exit 0
      ;;
    *)  break;;    # terminate while loop
  esac
  shift
done

if [ "$mach" = "pm-cpu" ]; then


    MODULE_FILE=$DIR/modules.pm-cpu.gnu
    export PETSC_DIR=/global/cfs/projectdirs/m4267/petsc/petsc_main/

    if [ "$config" -eq 1 ]; then
        export PETSC_ARCH=pm-cpu-hdf5_1_14_3-opt-64bit-gcc-11-2-0-fc2888174f5
    elif [ "$config" -eq 2 ]; then
        export PETSC_ARCH=pm-cpu-hdf5_1_14_3-debug-64bit-gcc-11-2-0-fc2888174f5
    elif [ "$config" -eq 3 ]; then
        export PETSC_ARCH=pm-cpu-opt-32bit-gcc-11-2-0-fc2888174f5
    fi

elif [ "$mach" = "pm-gpu" ]; then

    MODULE_FILE=$DIR/modules.pm-gpu.gnugpu
    export PETSC_DIR=/global/cfs/projectdirs/m4267/petsc/petsc_main/

    if [ "$config" -eq 1 ]; then
      export PETSC_ARCH=pm-gpu-hdf5_1_14_3-opt-64bit-gcc-11-2-0-fc2888174f5
    elif [ "$config" -eq 2 ]; then
      export PETSC_ARCH=pm-gpu-hdf5_1_14_3-debug-64bit-gcc-11-2-0-fc2888174f5
    elif [ "$config" -eq 3 ]; then
      export PETSC_ARCH=pm-gpu-opt-32bit-gcc-11-2-0-fc2888174f5
    else
      echo "Unsupported --config $config "
      display_help 
    fi

elif [ "$mach" = "frontier"  ]; then

  MODULE_FILE=$DIR/modules.frontier.gnugpu
  export PETSC_DIR=/lustre/orion/cli192/proj-shared/petsc

  if [ "$config" -eq 1 ]; then
      export PETSC_ARCH=frontier-gpu-opt-64bit-gcc-12-3-0-0d6defa7a01
  elif [ "$config" -eq 2 ]; then
      export PETSC_ARCH=frontier-gpu-debug-64bit-gcc-12-3-0-0d6defa7a01
  elif [ "$config" -eq 3 ]; then
      export PETSC_ARCH=frontier-gpu-opt-32bit-gcc-12-3-0-0d6defa7a01
  fi

#elif [ "$mach" = "aurora"  ]; then
#
#  MODULE_FILE=$DIR/modules.aurora.oneapi
#  export PETSC_DIR=/lus/gecko/projects/CSC250STMS07_CNDA/bishtgautam/petsc
#  if [ "$with64bit" -eq 0 ]; then
#    if [ "$with_debugging" -eq 0 ]; then
#      export PETSC_ARCH=aurora-opt-32bit-oneapi-ifx-fc288817
#    else
#       echo "On Aurora, --with-debugging 1 was selected, but PETSc has not been installed with debugging turned on."
#       exit 0
#    fi
#  else
#    if [ "$with_debugging" -eq 0 ]; then
#      export PETSC_ARCH=aurora-opt-64bit-oneapi-ifx-fc288817
#    else
#       echo "On Aurora, --with-debugging 1 was selected, but PETSc has not been installed with debugging turned on."
#       exit 0
#    fi
#  fi
#
else
  echo "Could not determine the machine. mach=$mach"
  display_help
  exit 0
fi

echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "Will source the following module file and set the following PETSc settings"
echo ""
echo "  source $MODULE_FILE"
echo "  export PETSC_DIR=$PETSC_DIR"
echo "  export PETSC_ARCH=$PETSC_ARCH"
echo ""
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
source $MODULE_FILE
export PETSC_DIR=$PETSC_DIR
export PETSC_ARCH=$PETSC_ARCH

