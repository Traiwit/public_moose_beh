
      subroutine DATE_AND_TIME(currentdate, currenttime)

      character*8, intent(out) :: currentdate
      character*10, intent(out) :: currenttime

      currentdate = "32.12.99"
      currenttime = "25:61:61.0"
      return
      end

      subroutine XPLB_EXIT
      write (*,*) "EXIT: XPLB_EXIT has been called."
      stop 1
      end
      
C     MPI related subroutines, data

      function GETCOMMUNICATOR()
      integer :: GETCOMMUNICATOR
      GETCOMMUNICATOR = 0
      return
      end function

      subroutine VGETRANK( KPROCESSNUM ) ! get process number
      integer, intent(out) :: KPROCESSNUM
      KPROCESSNUM = 1
      return
      end

      subroutine MPI_COMM_SIZE(mpi_comm, num_procs, ierr)
      integer, intent(in) :: mpi_comm
      integer, intent(out) :: num_procs, ierr
      num_procs = 0
      ierr = 0
      return
      end

      subroutine VGETNUMCPUS(num_procs)
      integer, intent(out) :: num_procs
      num_procs = 1
      return
      end

      subroutine MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
      dimension BUFFER(*)
      integer :: COUNT, DATATYPE, ROOT, COMM, IERROR
      return
      end

      subroutine MPI_SEND(
     &     data_to_send, send_count, send_type,
     &     destination_ID, tag, comm, ierr)
      dimension data_to_send(*)
      integer, intent(in) :: send_count, send_type
      integer, intent(in) :: destination_ID, tag, comm
      integer, intent(out) :: ierr
      ierr = 0
      return
      end

      subroutine MPI_RECV(
     &     received_data, receive_count, receive_type,
     &     sender_ID, tag, comm, status, ierr)
      dimension received_data(*)
      integer :: receive_count, receive_type
      integer :: sender_ID, tag, comm
      integer :: status, ierr
      ierr = 0
      return
      end

      subroutine MPI_TEST(REQUEST, FLAG, STATUS, IERROR)
      LOGICAL FLAG
      INTEGER REQUEST, STATUS(8), IERROR
C     INTEGER REQUEST, STATUS(MPI_STATUS_SIZE), IERROR; MPI_STATUS_SIZE = 8
      FLAG = .TRUE.
      return
      end

      subroutine MPI_IRECV(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, REQUEST, IERROR)
      DOUBLE PRECISION BUF(*)
      INTEGER COUNT, DATATYPE, SOURCE, TAG, COMM, REQUEST, IERROR
      IERROR = 0
      return
      end

      BLOCK DATA MPI_VALUES
      integer*4 MPI_STATUS_IGNORE
      common/mpi_f_status_ignore/MPI_STATUS_IGNORE
      DATA  MPI_STATUS_IGNORE / 1.0 /
      save /mpi_f_status_ignore/
      END

C     logfiles stuff....

      subroutine openLogFile(unit, filename)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: filename

      open(unit, FILE=trim(filename))
      return
      end


      subroutine openLogFileScratch(unit)
      integer, intent(in) :: unit
      open(unit,  STATUS='SCRATCH')
      return
      end


      subroutine closeLogFile(unit)
      integer, intent(in) :: unit
      close(unit)
      return
      end
