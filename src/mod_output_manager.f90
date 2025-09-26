module output_manager
  !=================================================================================
  ! Generic Data Output System
  ! Property-agnostic module for consistent data output in any format
  ! No knowledge of physics, chemistry, or specific property types
  !=================================================================================
  
  implicit none
  
  ! ========================= Output Format Constants ===========================
  
  integer, parameter :: OUTPUT_SCALAR = 1        ! Single number
  integer, parameter :: OUTPUT_VECTOR = 2        ! 1D array with optional coordinates
  integer, parameter :: OUTPUT_MATRIX = 3        ! 2D array
  integer, parameter :: OUTPUT_CUBE = 4          ! 3D regular grid (cube format)
  integer, parameter :: OUTPUT_POINTS = 5        ! Irregular 3D points
  integer, parameter :: OUTPUT_TABLE = 6         ! Multi-column tabular data
  
  ! ========================= File Format Constants ==============================
  
  character(len=*), parameter :: EXT_SCALAR = '.value'
  character(len=*), parameter :: EXT_VECTOR = '.dat'
  character(len=*), parameter :: EXT_MATRIX = '.matrix'
  character(len=*), parameter :: EXT_CUBE = '.cube'
  character(len=*), parameter :: EXT_POINTS = '.xyz'
  character(len=*), parameter :: EXT_TABLE = '.table'
  
  ! ========================= Error Codes =====================================
  
  integer, parameter :: STATUS_SUCCESS = 0
  integer, parameter :: STATUS_FILE_ERROR = 1
  integer, parameter :: STATUS_INVALID_FORMAT = 2
  integer, parameter :: STATUS_DATA_ERROR = 3
  integer, parameter :: STATUS_DIMENSION_ERROR = 4
  
  ! ========================= Data Structures =================================
  
  type :: output_config
    ! Basic identification
    character(len=128) :: title                 ! Free-form title for output
    character(len=256) :: filename              ! Output filename (auto-generated if empty)
    integer :: format_type                      ! OUTPUT_SCALAR, OUTPUT_VECTOR, etc.
    
    ! Formatting options
    integer :: precision                        ! Decimal places (default: 16)
    logical :: scientific_notation              ! Use exponential format
    logical :: include_header                   ! Include metadata header
    logical :: append_mode                      ! Append vs overwrite
    
    ! Grid parameters (for OUTPUT_CUBE only)
    real(8), dimension(3) :: origin             ! Grid origin [x,y,z]
    real(8), dimension(3) :: spacing            ! Grid spacing [dx,dy,dz]
    integer, dimension(3) :: dimensions         ! Grid size [nx,ny,nz]
    
    ! Column labels (for OUTPUT_TABLE, OUTPUT_VECTOR)
    character(len=32), allocatable :: column_labels(:)
    
    ! Metadata
    character(len=512) :: description           ! Free-form description
    character(len=256) :: units                 ! Physical units (optional)
    character(len=128) :: method                ! Computation method (optional)
    real(8) :: computation_time                 ! Timing info (optional)
    character(len=512) :: notes                 ! Additional notes
  end type output_config
  
  type :: data_container
    ! Scalar data
    real(8) :: scalar_value
    
    ! Vector data (1D)
    integer :: vector_size
    real(8), allocatable :: vector_data(:)
    real(8), allocatable :: vector_coords(:)    ! Optional x-coordinates
    
    ! Matrix data (2D)
    integer :: matrix_rows, matrix_cols
    real(8), allocatable :: matrix_data(:,:)
    
    ! Point cloud data (irregular 3D)
    integer :: n_points
    real(8), allocatable :: point_coords(:,:)   ! Shape: (3, n_points)
    real(8), allocatable :: point_values(:)     ! Values at each point
    
    ! Table data (multi-column)
    integer :: table_rows, table_cols
    real(8), allocatable :: table_data(:,:)     ! Shape: (table_rows, table_cols)
    
    ! Grid data (regular 3D - flattened)
    integer :: grid_size
    real(8), allocatable :: grid_values(:)      ! Flattened 3D data
  end type data_container
  
  ! ========================= Public Interface ================================
  
  public :: output_config, data_container
  public :: OUTPUT_SCALAR, OUTPUT_VECTOR, OUTPUT_MATRIX, OUTPUT_CUBE, OUTPUT_POINTS, OUTPUT_TABLE
  public :: STATUS_SUCCESS, STATUS_FILE_ERROR, STATUS_INVALID_FORMAT, STATUS_DATA_ERROR
  public :: create_output_config, write_data, cleanup_data
  
contains

  ! ========================= Configuration Factory ==========================
  
  function create_output_config(title, format_type, filename) result(config)
    character(len=*), intent(in) :: title
    integer, intent(in) :: format_type
    character(len=*), intent(in), optional :: filename
    type(output_config) :: config
    
    ! Initialize basic settings
    config%title = trim(title)
    config%format_type = format_type
    
    ! Set filename or auto-generate
    if (present(filename)) then
      config%filename = trim(filename)
    else
      config%filename = generate_filename(title, format_type)
    end if
    
    ! Set defaults
    config%precision = 16
    config%scientific_notation = .true.
    config%include_header = .true.
    config%append_mode = .false.
    
    ! Initialize grid parameters
    config%origin = [0.0d0, 0.0d0, 0.0d0]
    config%spacing = [1.0d0, 1.0d0, 1.0d0]
    config%dimensions = [1, 1, 1]
    
    ! Initialize metadata
    config%description = ''
    config%units = ''
    config%method = ''
    config%computation_time = 0.0d0
    config%notes = ''
    
  end function create_output_config
  
  ! ========================= Main Output Interface ===========================
  
  subroutine write_data(config, data, status)
    type(output_config), intent(in) :: config
    type(data_container), intent(in) :: data
    integer, intent(out), optional :: status
    
    integer :: unit_num, iostat, local_status
    logical :: file_exists
    
    local_status = STATUS_SUCCESS
    
    ! Validate inputs
    if (.not. validate_config(config)) then
      local_status = STATUS_INVALID_FORMAT
      if (present(status)) status = local_status
      return
    end if
    
    if (.not. validate_data(config, data)) then
      local_status = STATUS_DATA_ERROR
      if (present(status)) status = local_status
      return
    end if
    
    ! Handle file opening
    inquire(file=trim(config%filename), exist=file_exists)
    
    unit_num = get_free_unit()
    if (config%append_mode .and. file_exists) then
      open(unit=unit_num, file=trim(config%filename), status='old', &
           position='append', iostat=iostat)
    else
      open(unit=unit_num, file=trim(config%filename), status='replace', iostat=iostat)
    end if
    
    if (iostat /= 0) then
      local_status = STATUS_FILE_ERROR
      if (present(status)) status = local_status
      return
    end if
    
    ! Write header if requested
    if (config%include_header) then
      call write_header(unit_num, config)
    end if
    
    ! Write data based on format type
    select case(config%format_type)
      case(OUTPUT_SCALAR)
        call write_scalar_data(unit_num, config, data)
      case(OUTPUT_VECTOR)
        call write_vector_data(unit_num, config, data)
      case(OUTPUT_MATRIX)
        call write_matrix_data(unit_num, config, data)
      case(OUTPUT_CUBE)
        call write_cube_data(unit_num, config, data)
      case(OUTPUT_POINTS)
        call write_points_data(unit_num, config, data)
      case(OUTPUT_TABLE)
        call write_table_data(unit_num, config, data)
      case default
        local_status = STATUS_INVALID_FORMAT
    end select
    
    close(unit_num)
    
    if (present(status)) status = local_status
    
  end subroutine write_data
  
  ! ========================= Format-Specific Writers =========================
  
  subroutine write_scalar_data(unit_num, config, data)
    integer, intent(in) :: unit_num
    type(output_config), intent(in) :: config
    type(data_container), intent(in) :: data
    
    character(len=32) :: format_str
    
    ! Create format string
    call create_real_format(config, format_str)
    
    write(unit_num, '(A)') '# Scalar Value'
    write(unit_num, '(A)') '# ----------------------------------------'
    write(unit_num, format_str) data%scalar_value
    
  end subroutine write_scalar_data
  
  subroutine write_vector_data(unit_num, config, data)
    integer, intent(in) :: unit_num
    type(output_config), intent(in) :: config
    type(data_container), intent(in) :: data
    
    character(len=64) :: format_str
    integer :: i
    logical :: has_coords
    
    has_coords = allocated(data%vector_coords)
    
    ! Write column headers
    write(unit_num, '(A)') '# Vector Data'
    write(unit_num, '(A)') '# ----------------------------------------'
    
    if (allocated(config%column_labels)) then
      if (has_coords .and. size(config%column_labels) >= 2) then
        write(unit_num, '(A)') '# Columns: ' // trim(config%column_labels(1)) // '  ' // trim(config%column_labels(2))
      else if (.not. has_coords .and. size(config%column_labels) >= 1) then
        write(unit_num, '(A)') '# Column: ' // trim(config%column_labels(1))
      end if
    else
      if (has_coords) then
        write(unit_num, '(A)') '# Columns: X  Y'
      else
        write(unit_num, '(A)') '# Column: Value'
      end if
    end if
    
    write(unit_num, '(A)') '# ----------------------------------------'
    
    ! Create format string
    if (has_coords) then
      call create_real_format(config, format_str, 2)
    else
      call create_real_format(config, format_str, 1)
    end if
    
    ! Write data
    do i = 1, data%vector_size
      if (has_coords) then
        write(unit_num, format_str) data%vector_coords(i), data%vector_data(i)
      else
        write(unit_num, format_str) data%vector_data(i)
      end if
    end do
    
  end subroutine write_vector_data
  
  subroutine write_matrix_data(unit_num, config, data)
    integer, intent(in) :: unit_num
    type(output_config), intent(in) :: config
    type(data_container), intent(in) :: data
    
    character(len=128) :: format_str
    integer :: i, j
    
    write(unit_num, '(A)') '# Matrix Data'
    write(unit_num, '(A, I0, A, I0)') '# Dimensions: ', data%matrix_rows, ' x ', data%matrix_cols
    write(unit_num, '(A)') '# ----------------------------------------'
    
    ! Create format string for full row
    call create_real_format(config, format_str, data%matrix_cols)
    
    ! Write matrix data
    do i = 1, data%matrix_rows
      write(unit_num, format_str) (data%matrix_data(i, j), j = 1, data%matrix_cols)
    end do
    
  end subroutine write_matrix_data
  
  subroutine write_cube_data(unit_num, config, data)
    integer, intent(in) :: unit_num
    type(output_config), intent(in) :: config
    type(data_container), intent(in) :: data
    
    integer :: nx, ny, nz, i, j, k, idx
    
    nx = config%dimensions(1)
    ny = config%dimensions(2)
    nz = config%dimensions(3)
    
    ! Standard cube file format header
    write(unit_num, '(A)') trim(config%title)
    write(unit_num, '(A)') trim(config%description)
    
    ! Number of atoms and origin (dummy atom for pure data)
    write(unit_num, '(I5, 3F12.6)') 1, config%origin
    
    ! Voxel vectors
    write(unit_num, '(I5, 3F12.6)') nx, config%spacing(1), 0.0d0, 0.0d0
    write(unit_num, '(I5, 3F12.6)') ny, 0.0d0, config%spacing(2), 0.0d0
    write(unit_num, '(I5, 3F12.6)') nz, 0.0d0, 0.0d0, config%spacing(3)
    
    ! Dummy atom at origin
    write(unit_num, '(I5, 4F12.6)') 1, 0.0d0, 0.0d0, 0.0d0, 0.0d0
    
    ! Data values (x outermost, z innermost)
    idx = 0
    do i = 1, nx
      do j = 1, ny
        do k = 1, nz
          idx = idx + 1
          if (idx <= data%grid_size) then
            write(unit_num, '(ES13.5)', advance='no') data%grid_values(idx)
          else
            write(unit_num, '(ES13.5)', advance='no') 0.0d0
          end if
          
          ! Line breaks every 6 values
          if (mod(k, 6) == 0 .or. k == nz) then
            write(unit_num, *)
          end if
        end do
      end do
    end do
    
  end subroutine write_cube_data
  
  subroutine write_points_data(unit_num, config, data)
    integer, intent(in) :: unit_num
    type(output_config), intent(in) :: config
    type(data_container), intent(in) :: data
    
    character(len=64) :: format_str
    integer :: i
    
    write(unit_num, '(A)') '# Point Cloud Data'
    write(unit_num, '(A, I0)') '# Number of points: ', data%n_points
    write(unit_num, '(A)') '# ----------------------------------------'
    
    if (allocated(config%column_labels) .and. size(config%column_labels) >= 4) then
      write(unit_num, '(A)') '# Columns: ' // trim(config%column_labels(1)) // '  ' // &
                             trim(config%column_labels(2)) // '  ' // &
                             trim(config%column_labels(3)) // '  ' // &
                             trim(config%column_labels(4))
    else
      write(unit_num, '(A)') '# Columns: X  Y  Z  Value'
    end if
    
    write(unit_num, '(A)') '# ----------------------------------------'
    
    ! Create format string for 4 columns
    call create_real_format(config, format_str, 4)
    
    ! Write point data
    do i = 1, data%n_points
      write(unit_num, format_str) data%point_coords(1, i), data%point_coords(2, i), &
                                  data%point_coords(3, i), data%point_values(i)
    end do
    
  end subroutine write_points_data
  
  subroutine write_table_data(unit_num, config, data)
    integer, intent(in) :: unit_num
    type(output_config), intent(in) :: config
    type(data_container), intent(in) :: data
    
    character(len=256) :: format_str, header_line
    integer :: i, j
    
    write(unit_num, '(A)') '# Tabular Data'
    write(unit_num, '(A, I0, A, I0)') '# Dimensions: ', data%table_rows, ' rows x ', data%table_cols, ' columns'
    write(unit_num, '(A)') '# ----------------------------------------'
    
    ! Write column headers if available
    if (allocated(config%column_labels) .and. size(config%column_labels) >= data%table_cols) then
      header_line = '# Columns:'
      do j = 1, data%table_cols
        header_line = trim(header_line) // '  ' // trim(config%column_labels(j))
      end do
      write(unit_num, '(A)') trim(header_line)
      write(unit_num, '(A)') '# ----------------------------------------'
    end if
    
    ! Create format string for all columns
    call create_real_format(config, format_str, data%table_cols)
    
    ! Write table data
    do i = 1, data%table_rows
      write(unit_num, format_str) (data%table_data(i, j), j = 1, data%table_cols)
    end do
    
  end subroutine write_table_data
  
  ! ========================= Utility Functions ===============================
  
  subroutine write_header(unit_num, config)
    integer, intent(in) :: unit_num
    type(output_config), intent(in) :: config
    
    character(len=24) :: timestamp
    
    ! Get timestamp
    call date_and_time(date=timestamp(1:8), time=timestamp(9:18))
    
    write(unit_num, '(A)') '# =========================================='
    write(unit_num, '(A)') '# ' // trim(config%title)
    write(unit_num, '(A)') '# Generated: ' // timestamp(1:4) // '-' // timestamp(5:6) // '-' // timestamp(7:8) // &
                          ' ' // timestamp(9:10) // ':' // timestamp(11:12) // ':' // timestamp(13:14)
    write(unit_num, '(A)') '# =========================================='
    
    if (len_trim(config%description) > 0) then
      write(unit_num, '(A)') '# Description: ' // trim(config%description)
    end if
    
    if (len_trim(config%units) > 0) then
      write(unit_num, '(A)') '# Units: ' // trim(config%units)
    end if
    
    if (len_trim(config%method) > 0) then
      write(unit_num, '(A)') '# Method: ' // trim(config%method)
    end if
    
    if (config%computation_time > 0.0d0) then
      write(unit_num, '(A, F12.6, A)') '# Computation time: ', config%computation_time, ' seconds'
    end if
    
    if (len_trim(config%notes) > 0) then
      write(unit_num, '(A)') '# Notes: ' // trim(config%notes)
    end if
    
  end subroutine write_header
  
  subroutine create_real_format(config, format_str, n_columns)
    type(output_config), intent(in) :: config
    character(len=*), intent(out) :: format_str
    integer, intent(in), optional :: n_columns
    
    integer :: cols, field_width
    character(len=8) :: exp_format
    
    cols = 1
    if (present(n_columns)) cols = n_columns
    
    field_width = config%precision + 8  ! Extra space for sign, decimal, exponent
    
    if (config%scientific_notation) then
      write(exp_format, '(A,I0,A)') 'E', config%precision, 'E3'
      write(format_str, '(A,I0,A,I0,A)') '(', cols, 'ES', field_width, '.', exp_format, ')'
    else
      write(format_str, '(A,I0,A,I0,A,I0,A)') '(', cols, 'F', field_width, '.', config%precision, ')'
    end if
    
  end subroutine create_real_format
  
  function generate_filename(title, format_type) result(filename)
    character(len=*), intent(in) :: title
    integer, intent(in) :: format_type
    character(len=256) :: filename
    
    character(len=128) :: clean_title
    character(len=8) :: extension
    
    ! Clean the title for use as filename
    clean_title = sanitize_filename(title)
    
    ! Choose extension based on format
    select case(format_type)
      case(OUTPUT_SCALAR)
        extension = EXT_SCALAR
      case(OUTPUT_VECTOR)
        extension = EXT_VECTOR
      case(OUTPUT_MATRIX)
        extension = EXT_MATRIX
      case(OUTPUT_CUBE)
        extension = EXT_CUBE
      case(OUTPUT_POINTS)
        extension = EXT_POINTS
      case(OUTPUT_TABLE)
        extension = EXT_TABLE
      case default
        extension = '.dat'
    end select
    
    filename = trim(clean_title) // trim(extension)
    
  end function generate_filename
  
  function sanitize_filename(name) result(clean_name)
    character(len=*), intent(in) :: name
    character(len=128) :: clean_name
    integer :: i
    character :: c
    
    clean_name = ''
    
    do i = 1, min(len_trim(name), 128)
      c = name(i:i)
      
      ! Replace problematic characters with underscores
      if (c == ' ' .or. c == '/' .or. c == '\' .or. c == ':' .or. &
          c == '*' .or. c == '?' .or. c == '"' .or. c == '<' .or. &
          c == '>' .or. c == '|' .or. c == '(' .or. c == ')' .or. &
          c == '[' .or. c == ']' .or. c == '{' .or. c == '}') then
        clean_name = trim(clean_name) // '_'
      else if ((c >= 'a' .and. c <= 'z') .or. (c >= 'A' .and. c <= 'Z') .or. &
               (c >= '0' .and. c <= '9') .or. c == '_' .or. c == '-' .or. c == '.') then
        clean_name = trim(clean_name) // c
      end if
    end do
    
    ! Ensure name is not empty
    if (len_trim(clean_name) == 0) then
      clean_name = 'output'
    end if
    
  end function sanitize_filename
  
  function get_free_unit() result(unit_num)
    integer :: unit_num
    logical :: is_open
    
    do unit_num = 10, 99
      inquire(unit=unit_num, opened=is_open)
      if (.not. is_open) exit
    end do
    
  end function get_free_unit
  
  ! ========================= Validation Functions ============================
  
  function validate_config(config) result(is_valid)
    type(output_config), intent(in) :: config
    logical :: is_valid
    
    is_valid = .true.
    
    ! Check format type
    if (config%format_type < 1 .or. config%format_type > 6) then
      is_valid = .false.
      return
    end if
    
    ! Check precision
    if (config%precision < 1 .or. config%precision > 20) then
      is_valid = .false.
      return
    end if
    
    ! Check grid dimensions for cube format
    if (config%format_type == OUTPUT_CUBE) then
      if (any(config%dimensions <= 0) .or. any(config%spacing <= 0.0d0)) then
        is_valid = .false.
        return
      end if
    end if
    
  end function validate_config
  
  function validate_data(config, data) result(is_valid)
    type(output_config), intent(in) :: config
    type(data_container), intent(in) :: data
    logical :: is_valid
    
    is_valid = .true.
    
    select case(config%format_type)
      case(OUTPUT_SCALAR)
        ! No special validation needed for scalar
        
      case(OUTPUT_VECTOR)
        if (.not. allocated(data%vector_data) .or. data%vector_size <= 0) then
          is_valid = .false.
          return
        end if
        if (allocated(data%vector_coords)) then
          if (size(data%vector_coords) /= data%vector_size) then
            is_valid = .false.
            return
          end if
        end if
        
      case(OUTPUT_MATRIX)
        if (.not. allocated(data%matrix_data) .or. data%matrix_rows <= 0 .or. data%matrix_cols <= 0) then
          is_valid = .false.
          return
        end if
        
      case(OUTPUT_CUBE)
        if (.not. allocated(data%grid_values) .or. data%grid_size <= 0) then
          is_valid = .false.
          return
        end if
        
      case(OUTPUT_POINTS)
        if (.not. allocated(data%point_coords) .or. .not. allocated(data%point_values) .or. data%n_points <= 0) then
          is_valid = .false.
          return
        end if
        
      case(OUTPUT_TABLE)
        if (.not. allocated(data%table_data) .or. data%table_rows <= 0 .or. data%table_cols <= 0) then
          is_valid = .false.
          return
        end if
        
    end select
    
  end function validate_data
  
  ! ========================= Memory Management ===============================
  
  subroutine cleanup_data(data)
    type(data_container), intent(inout) :: data
    
    if (allocated(data%vector_data)) deallocate(data%vector_data)
    if (allocated(data%vector_coords)) deallocate(data%vector_coords)
    if (allocated(data%matrix_data)) deallocate(data%matrix_data)
    if (allocated(data%point_coords)) deallocate(data%point_coords)
    if (allocated(data%point_values)) deallocate(data%point_values)
    if (allocated(data%table_data)) deallocate(data%table_data)
    if (allocated(data%grid_values)) deallocate(data%grid_values)
    
  end subroutine cleanup_data

end module output_manager
